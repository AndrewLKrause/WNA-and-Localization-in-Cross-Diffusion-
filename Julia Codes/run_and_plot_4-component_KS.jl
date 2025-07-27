t_start = time()

include("cross_diffusion_1D.jl")
include("plot_solution.jl")

N = 2
L = 20.0
Nx = 200
tspan = (0.0, 100.0)
Ncomp = 4

c  =  5.0

a  =  0.1
r_1  =  1.0
b  =  0.1
r_2 = 1.0
r  =  0.1

# Nonlinear diffusion matrix function (in-place)
function D_func!(D, u)
    D[1,1] = 1.0
    D[1,2] = -c*u[1] ./(1+u[1] .*u[1])
    D[1,3] = 0.0
    D[1,4] = 0.0
    D[2,1] = 0.0
    D[2,2] = 1.0
    D[2,3] = 0.0
    D[2,4] = 0.0
    D[3,1] = 0.0
    D[3,2] = 0.0
    D[3,3] = 1.0
    D[3,4] = 0.0
    D[4,1] = 0.0
    D[4,2] = 0.0
    D[4,3] = 0.0
    D[4,4] = 1.0
end

# Reaction function (in-place)
function R_func!(R, u)
    R[1] = r*u[1]*(1-u[1])
    R[2] = u[1]-a*u[2]-r_2*u[2]*u[3]+b*u[4]
    R[3] = -r_2*u[2]*u[3]+(r_1+b)*u[4]+u[1]*(1-u[3])
    R[4] = r_2*u[2]*u[3]-(r_1+b)*u[4]
end

# Initial condition function
function u0_func(x)
    [1+0.1*cos(2*pi*x/L),1,1,1]
end


x, sol = cross_diffusion_1D(D_func!, R_func!, u0_func, Ncomp, L, Nx, tspan)

plot_solution(x, sol; component=1, times=[0.0, 10.0, 50.0, 100.0])

elapsed = time() - t_start
println("Elapsed time: $elapsed seconds")