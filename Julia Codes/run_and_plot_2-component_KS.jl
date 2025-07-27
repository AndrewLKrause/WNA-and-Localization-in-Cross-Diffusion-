t_start = time()


N = 2
L = 10.0
Nx = 100
tspan = (0.0, 100.0)

r = 1.0
K = 2.0
b = 1.0
a = 0.5
chi  =  0.589

# Nonlinear diffusion matrix function (in-place)
function D_func!(D, u)
    D[1,1] = 0.1
    D[1,2] = -chi*u[1]
    D[2,1] = 0.0
    D[2,2] = 0.1
end

# Reaction function (in-place)
function R_func!(R, u)
    R[1] = r*u[1]*(1-u[1]/K)
    R[2] = a*u[1]-b*u[2]
end

# Initial condition function
function u0_func(x)
    [K+0.1*cos(2*pi*x/L), K*a/b]
end


x, sol = cross_diffusion_1D(D_func!, R_func!, u0_func, Ncomp, L, Nx, tspan)

plot_solution(x, sol; component=1, times=[0.0, 10.0, 50.0, 100.0])

elapsed = time() - t_start
println("Elapsed time: $elapsed seconds")