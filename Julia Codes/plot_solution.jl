using Plots

function plot_solution(x, sol; component=1, times)
    Nx = length(x)
    N = length(sol.u[1]) ÷ Nx  # number of components

    plt = plot(title="Component $component over time", xlabel="x", ylabel="u[$component]")
    for t in times
        u_t = sol(t)
        u_comp = [u_t[(i-1)*N + component] for i in 1:Nx]
        plot!(x, u_comp, label="t = $t")
    end
    display(plt)
end