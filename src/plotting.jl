function plot_results(time_array::AbstractArray, results::AbstractArray, parameters::AbstractArray, type_of_results::AbstractString)
    """
    General function to save and plot the non-projected results.
    """
    ω = parameters[1] 
    plot_position(time_array, real.(results[:,2]),"Real", type_of_results)
    plot_position(time_array, imag.(results[:,2]),"Imag", type_of_results)
    plot_current(time_array, real.(results[:,1]), "Real", type_of_results)
    plot_current(time_array, imag.(results[:,1]), "Imag", type_of_results)
    plot_fourier(time_array, results[:,1], ω, type_of_results)
    plot_norm(time_array, norm.(results[:,3]) .^ 2, type_of_results)
    plot_groundstate(time_array, real.(results[:,4]), type_of_results)
    return nothing
end
    
function plot_eigen(eigenvalues, gap::T, E_exciton::T, flag::AbstractString, type_of_calc::AbstractString) where T <: Number
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Eigenvalues", xlabel = "index of eigs", ylabel = "Energy [au]")
    indexes = [i_ for i_ in eachindex(eigenvalues)]
    scatter!(ax,indexes, eigenvalues, label = "Eigenvalues")
    axislegend()
    save("$current_path/plots/eigenvalues_$(flag)_$(type_of_calc).pdf", fig)
    save("$current_path/plots/eigenvalues_$(flag)_$(type_of_calc).png", fig)
    outfile = "$current_path/results/eigenvalues_$(flag)_$(type_of_calc).dat"
    open(outfile, "w") do f
        for (index,eig) in enumerate(eigenvalues)
          println(f, index, "\t", eig)
        end
    end

    if flag == "Real"
        fig_2 = Figure()
        ax_2 = Axis(fig_2[1,1], title = "Eigenvalues", ylabel = "Energy [au]",
            xticklabelsvisible = false, xticksvisible = false)
        hlines!(ax_2, eigenvalues, xmax = 0.75, xmin = 0.25, label = "Eigenvalues", color = :firebrick)
        hlines!(ax_2, -gap, label = L"$-\Delta_{hBN}$", color = :gray12, linestyle = :dash)
        hlines!(ax_2, -E_exciton, label = L"$-E_{g}$", color = :lightblue4, linestyle = :dashdot)
        ylims!(minimum(eigenvalues) - 0.01, gap)
        axislegend()
        save("$current_path/plots/eigenvalues_small_$(type_of_calc).pdf", fig_2)
        save("$current_path/plots/eigenvalues_small_$(type_of_calc).png", fig_2)
    end
    return nothing
end 

function plot_potential(xs::AbstractArray, gap::T, E_exciton::T, a_rydbergs::T, η_rydbergs::T, a_no_rydbergs::T, η_no_rydbergs::T) where T <: Number
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Potential energy", xlabel = L"$x$ [au]", ylabel = L"$V(x)$ [au]")
    lines!(ax, xs, Vatom.(xs, a_rydbergs, η_rydbergs), label = "Rydbergs")
    lines!(ax, xs, Vatom.(xs, a_no_rydbergs, η_no_rydbergs), label = "No Rydbergs")
    hlines!(ax, -eV_to_au(gap), xmin = minimum(xs), xmax = maximum(xs), label = L"$\Delta_{hBN}$", color = :grey)
    hlines!(ax, -eV_to_au(E_exciton), xmin = minimum(xs), xmax = maximum(xs), label = L"$\Delta_{hBN}$", color = :grey, linestyle = :dash)
    axislegend()
    save("$current_path/plots/potential.pdf", fig)
    save("$current_path/plots/potential.png", fig)
    outfile = "$current_path/results/potential_parameters.dat"
    labels = ["a_rydbergs", "η_rydbergs", "a_no_rydbergs", "η_no_rydbergs"]
    values = [a_rydbergs, η_rydbergs, a_no_rydbergs, η_no_rydbergs]
    open(outfile, "w") do f
        for (label, value) in zip(labels, values)
          println(f, label, "\t", value)
        end
    end

    return nothing
end

function plot_electric_field(ts::AbstractArray, E_t::AbstractArray)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Electric field", xlabel = L"$t$ [au]", ylabel = L"$E(t)$ [au]")
    lines!(ax, ts, E_t, label = "Electric field")
    axislegend()
    save("$current_path/plots/electric_field.pdf", fig)
    save("$current_path/plots/electric_field.png", fig)
    return nothing
end

function plot_current(ts::AbstractArray, current::AbstractArray, flag::AbstractString, type_of_calc::AbstractString)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "$flag part of current", xlabel = L"$t$ [au]", ylabel = L"$j(t)$ [au]")
    lines!(ax, ts, current, label = L"$j(t)$")
    axislegend()
    save("$current_path/plots/current_$(flag)_$(type_of_calc).pdf", fig)
    save("$current_path/plots/current_$(flag)_$(type_of_calc).png", fig)  
    outfile = "$current_path/results/current_$(flag)_$(type_of_calc).dat"
    open(outfile, "w") do f
        for (t, jt) in zip(ts, current)
          println(f, t, "\t", jt)
        end
    end
    return nothing
end

function plot_position(ts::AbstractArray, position::AbstractArray, flag::AbstractString, type_of_calc::AbstractString)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "$flag part of position", xlabel = L"$t$ [au]", ylabel = L"$r(t)$ [au]")
    lines!(ax, ts, position, label = L"$r(t)$")
    axislegend()
    save("$current_path/plots/position_$(flag)_$(type_of_calc).pdf", fig)
    save("$current_path/plots/position_$(flag)_$(type_of_calc).png", fig)  
    outfile = "$current_path/results/position_$(flag)_$(type_of_calc).dat"
    open(outfile, "w") do f
        for (t, rt) in zip(ts, position)
          println(f, t, "\t", rt)
        end
    end
    return nothing
end


function plot_norm(ts::AbstractArray, norm::AbstractArray, type_of_calc::AbstractString)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Norm", xlabel = L"$t$ [au]", ylabel = L"$\langle \Psi(t)| \Psi(t)\rangle $")
    lines!(ax, ts, norm, label = L"$\langle \Psi(t) | \Psi(t) \rangle $")
    axislegend()
    save("$current_path/plots/norm_$(type_of_calc).pdf", fig)
    save("$current_path/plots/norm_$(type_of_calc).png", fig)  
    outfile = "$current_path/results/norm_$(type_of_calc).dat"
    open(outfile, "w") do f
        for (t, norm_t) in zip(ts, norm)
          println(f, t, "\t", norm_t)
        end
    end
    return nothing
end

function plot_groundstate(ts::AbstractArray, norm::AbstractArray, type_of_calc::AbstractString)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Norm", xlabel = L"$t$ [au]", ylabel = L"$\langle \Psi_0| \Psi(t)\rangle $")
    lines!(ax, ts, norm, label = L"$\langle \Psi_0 | \Psi(t) \rangle $")
    axislegend()
    save("$current_path/plots/ground_state_$(type_of_calc).pdf", fig)
    save("$current_path/plots/ground_state_$(type_of_calc).png", fig)  
    outfile = "$current_path/results/ground_state_$(type_of_calc).dat"
    open(outfile, "w") do f
        for (t, norm_t) in zip(ts, norm)
          println(f, t, "\t", norm_t)
        end
    end
    return nothing
end


function plot_fourier(time_array::AbstractArray, current::AbstractArray, ω::Real, type_of_calc::AbstractString)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], yscale = log10, title = "Spectrum", xlabel = "Harmonic Number", ylabel = "Spectrum [a.u.]")
    ω_array, spectrum = harmonic_analysis(time_array, real.(current), ω)
    ω_array ./= ω
    lines!(ax, ω_array, spectrum, label = "Spectrum")
    ylims!(ax, [1e-16, 1e-2])
    xlims!(ax, 0, 100)
    axislegend()
    save("$current_path/plots/spectrum_$(type_of_calc).pdf", fig)
    save("$current_path/plots/spectrum_$(type_of_calc).png", fig)  
    outfile = "$current_path/results/spectrum_$(type_of_calc).dat"
    open(outfile, "w") do f
        for (freq, spec) in zip(ω_array, spectrum)
          println(f, freq, "\t", spec)
        end
    end
    return nothing
end


###################################################################################################
#############################  OLD VERSIONS OF PROJECTIONS OF PLOTTING ############################
###################################################################################################

#= function plot_results_rydbergs(time_array::AbstractArray, results_rydbergs::AbstractArray, parameters::AbstractArray)
    """
    General function to save and plot the projected results.
    """
    ω = parameters[1] 
    plot_position_rydbergs(time_array, real.(results_rydbergs[:,2]),"Real")
    plot_position_rydbergs(time_array, imag.(results_rydbergs[:,2]),"Imag")
    plot_current_rydbergs(time_array, real.(results_rydbergs[:,1]), "Real")
    plot_current_rydbergs(time_array, imag.(results_rydbergs[:,1]), "Imag")
    plot_fourier_rydbergs(time_array, results_rydbergs[:,1], ω)
    plot_norm_rydbergs(time_array, norm.(results_rydbergs[:,3]) .^ 2)
    plot_groundstate_rydbergs(time_array, real.(results_rydbergs[:,4]))
    return nothing
end

function plot_eigen(eigenvalues, gap::T, E_exciton::T) where T <: Number
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Eigenvalues", xlabel = "index of eigs", ylabel = "Energy [au]")
    indexes = [i_ for i_ in eachindex(eigenvalues)]
    scatter!(ax,indexes, eigenvalues, label = "Eigenvalues")
    axislegend()
    save("$current_path/plots/eigenvalues.pdf", fig)
    save("$current_path/plots/eigenvalues.png", fig)
    outfile = "$current_path/results/eigenvalues.dat"
    open(outfile, "w") do f
        for (index,eig) in enumerate(eigenvalues)
          println(f, index, "\t", eig)
        end
    end

    fig_2 = Figure()
    ax_2 = Axis(fig_2[1,1], title = "Eigenvalues", ylabel = "Energy [au]",
        xticklabelsvisible = false, xticksvisible = false)
    hlines!(ax_2, eigenvalues, xmax = 0.75, xmin = 0.25, label = "Eigenvalues", color = :firebrick)
    hlines!(ax_2, -gap, label = L"$-\Delta_{hBN}$", color = :gray12, linestyle = :dash)
    hlines!(ax_2, -E_exciton, label = L"$-E_{g}$", color = :lightblue4, linestyle = :dashdot)
    ylims!(minimum(eigenvalues) - 0.01, gap)
    axislegend()
    save("$current_path/plots/eigenvalues_small.pdf", fig_2)
    save("$current_path/plots/eigenvalues_small.png", fig_2)
    return nothing
end 

function plot_fourier_rydbergs(time_array::AbstractArray, current::AbstractArray, ω::Real)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], yscale = log10, title = "Spectrum", xlabel = "Harmonic Number", ylabel = "Spectrum [a.u.]")
    ω_array, spectrum = harmonic_analysis(time_array, real.(current), ω)
    ω_array ./= ω
    lines!(ax, ω_array, spectrum, label = "Spectrum")
    ylims!(ax, [1e-16, 1e-7])
    xlims!(ax, 0, 50)
    axislegend()
    save("$current_path/plots/spectrum_rydbergs.pdf", fig)
    save("$current_path/plots/spectrum_rydbergs.png", fig)  
    outfile = "$current_path/results/spectrum_rydbergs.dat"
    open(outfile, "w") do f
        for (freq, spec) in zip(ω_array, spectrum)
          println(f, freq, "\t", spec)
        end
    end
    return nothing
end

function plot_groundstate_rydbergs(ts::AbstractArray, norm::AbstractArray)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Norm", xlabel = L"$t$ [au]", ylabel = L"$\langle \Psi_0| \Psi(t)\rangle $")
    lines!(ax, ts, Float16.(norm), label = L"$\langle \Psi_0 | \Psi(t) \rangle $")
    axislegend()
    save("$current_path/plots/ground_state_rydbergs.pdf", fig)
    save("$current_path/plots/ground_state_rydbergs.png", fig)  
    outfile = "$current_path/results/ground_state_rydbergs.dat"
    open(outfile, "w") do f
        for (t, norm_t) in zip(ts, norm)
          println(f, t, "\t", norm_t)
        end
    end
    return nothing
end

function plot_norm_rydbergs(ts::AbstractArray, norm::AbstractArray)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "Norm", xlabel = L"$t$ [au]", ylabel = L"$\langle \Psi(t)|\Psi(t)\rangle$")
    lines!(ax, ts, Float16.(norm), label = L"$\langle \Psi(t)|\Psi(t)\rangle$")
    axislegend()
    save("$current_path/plots/norm_rydbergs.pdf", fig)
    save("$current_path/plots/norm_rydbergs.png", fig)  
    outfile = "$current_path/results/norm_rydbergs.dat"
    open(outfile, "w") do f
        for (t, norm_t) in zip(ts, norm)
          println(f, t, "\t", norm_t)
        end
    end
    return nothing
end

function plot_position_rydbergs(ts::AbstractArray, position::AbstractArray, flag::AbstractString)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "$flag part of position", xlabel = L"$t$ [au]", ylabel = L"$r(t)$ [au]")
    lines!(ax, ts, position, label = L"$r(t)$")
    axislegend()
    save("$current_path/plots/position_$(flag)_rydbergs.pdf", fig)
    save("$current_path/plots/position_$(flag)_rydbergs.png", fig)  
    outfile = "$current_path/results/position_$(flag)_rydbergs.dat"
    open(outfile, "w") do f
        for (t, rt) in zip(ts, position)
          println(f, t, "\t", rt)
        end
    end
    return nothing
end
function plot_current_rydbergs(ts::AbstractArray, current::AbstractArray, flag::AbstractString)
    current_path = pwd()
    fig = Figure()
    ax = Axis(fig[1,1], title = "$flag part of current", xlabel = L"$t$ [au]", ylabel = L"$j(t)$ [au]")
    lines!(ax, ts, current, label = L"$j(t)$")
    #ylims!(ax, [-0.5e20, 0.5e20])
    #xlims!(ax, [minimum(ts), maximum(ts)])
    axislegend()
    save("$current_path/plots/current_$(flag)_rydbergs.pdf", fig)
    save("$current_path/plots/current_$(flag)_rydbergs.png", fig)
    outfile = "$current_path/results/current_$(flag)_rydbergs.dat"
    open(outfile, "w") do f
        for (t, jt) in zip(ts, current)
          println(f, t, "\t", jt)
        end
    end
    return nothing
end
 =#