using DifferentialEquations
using LinearAlgebra
using Plots
using FiniteDiff

# Main PDE solver
function cross_diffusion_1D(D_func!, R_func!, u0_func, Ncomp, L, Nx, tspan)
    dx = L / (Nx - 1)
    x = range(0, stop=L, length=Nx)
    u0 = zeros(Ncomp * Nx)
    for i in 1:Nx
        u0[(i-1)*Ncomp+1:i*Ncomp] .= u0_func(x[i])
    end

    function dudt!(du, u, p, t)
        fill!(du, 0.0)
        D = zeros(Ncomp, Ncomp)
        R = zeros(Ncomp)
        uL = zeros(Ncomp)
        uC = zeros(Ncomp)
        uR = zeros(Ncomp)
        gradL = zeros(Ncomp)
        gradR = zeros(Ncomp)
        fluxL = zeros(Ncomp)
        fluxR = zeros(Ncomp)

        for i in 2:Nx-1
            @views uL .= u[(i-2)*Ncomp+1:(i-1)*Ncomp]
            @views uC .= u[(i-1)*Ncomp+1:i*Ncomp]
            @views uR .= u[i*Ncomp+1:(i+1)*Ncomp]

            @. gradL = (uC - uL) / dx
            @. gradR = (uR - uC) / dx

            D_func!(D, 0.5 .* (uL .+ uC)); mul!(fluxL, D, gradL)
            D_func!(D, 0.5 .* (uC .+ uR)); mul!(fluxR, D, gradR)

            R_func!(R, uC)
            @views du[(i-1)*Ncomp+1:i*Ncomp] .= (fluxR - fluxL) / dx .+ R
        end

        # Left boundary (Neumann)
        @views uC .= u[1:Ncomp]
        @views uR .= u[Ncomp+1:2*Ncomp]
        @. gradR = (uR - uC) / dx
        D_func!(D, 0.5 .* (uC .+ uR)); mul!(fluxR, D, gradR)
        fill!(fluxL, 0.0)
        R_func!(R, uC)
        @views du[1:Ncomp] .= (fluxR - fluxL) / dx .+ R

        # Right boundary (Neumann)
        @views uL .= u[(Nx-2)*Ncomp+1:(Nx-1)*Ncomp]
        @views uC .= u[(Nx-1)*Ncomp+1:end]
        @. gradL = (uC - uL) / dx
        D_func!(D, 0.5 .* (uL .+ uC)); mul!(fluxL, D, gradL)
        fill!(fluxR, 0.0)
        R_func!(R, uC)
        @views du[(Nx-1)*Ncomp+1:end] .= (fluxR - fluxL) / dx .+ R
    end

    prob = ODEProblem(dudt!, u0, tspan)
    sol = solve(prob, Rodas5(autodiff=false); reltol=1e-6, abstol=1e-6)
    return x, sol
end