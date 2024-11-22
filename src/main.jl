function setting_up_folders!()
    """To create all the folders where the data and the plots are being save"""
    current_path = pwd()
    rm("$current_path/plots", recursive = true, force = true)
    rm("$current_path/results", recursive = true, force = true)
    mkpath("$current_path/plots")
    mkpath("$current_path/results")
    return nothing
end

function solve_dynamics(input_file::AbstractString)
    """
    Main function of the module. It solves the whole dynamics while outputting observables and figures
    It only needs a dictionary filled with the parameters we want.    
    """
    println("="^50)
    println("Starting to solve the dynamics of the 1d Hydrogen atom!!!")
    # Set up of the folders
    setting_up_folders!()

    println(" · Reading input...")
    input = ReadInput(input_file)

    laser_type = input.model["laser type"]
    potential_type = input.model["potential type"]

    ω = input.physical_params["ω"]
    E0 = input.physical_params["E0"]
    τ = input.physical_params["τ"]
    gap_hBN = input.physical_params["hBN_gap"]
    E_exciton = input.physical_params["E_exciton"]

    projected_bool = input.calculation_details["projected_bool"]
    timeprop_bool = input.calculation_details["timeprop_bool"]
    optim_bool_no_ryd = input.calculation_details["optim_bool_no_ryd"]
    optim_bool_ryd = input.calculation_details["optim_bool_ryd"]
    php_calc_bool = input.calculation_details["php_calc_bool"]

    dx = input.numerical_params["dx"]
    xmax = input.numerical_params["xmax"]
    dt = input.numerical_params["time_spacing"]
    # Parameters from previous calculations...
    a_no_rydbergs_prev = input.numerical_params["a_no_rydbergs"]
    η_no_rydbergs_prev = input.numerical_params["η_no_rydbergs"]
    a_rydbergs_prev = input.numerical_params["a_rydbergs"]
    η_rydbergs_prev = input.numerical_params["η_rydbergs"]
    flush(stdout)

    print("---> Input read! "); printstyled(" ✓\n", color = :light_green)
    println(" · Constructing the time-independent quantities...")
    x_array = collect(-xmax:dx:xmax) 
    println(" · The dimension of the Hillbert space used is ",length(x_array))
    if optim_bool_no_ryd == true
        H0, a_no_rydbergs, η_no_rydbergs = hamiltonian_optimization(x_array, -gap_hBN, 0.0, a_no_rydbergs_prev, η_no_rydbergs_prev)
    else
        println(" · We are not optimizing!!!")
        H0 = hamiltonian_no_optimization(x_array, a_no_rydbergs_prev, η_no_rydbergs_prev)
        a_no_rydbergs = a_no_rydbergs_prev
        η_no_rydbergs = η_no_rydbergs_prev
    end 
    r_operator = R_operator(x_array)
    print("--> Hamiltonian and R operator built! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)
    println(" · Getting the spectrum of H ...")
    flush(stdout)
    #eigenvalues, eigenvectors = eigs(H0, nev = 1000, which = :SR)
    eigenvalues, eigenvectors = eigen(Array(H0))
    print("--> Spectrum obtained! "); printstyled(" ✓\n", color = :light_green)
    println(" · Plotting the spectrum of H ...")
    flush(stdout)
    plot_eigen(real.(eigenvalues), gap_hBN, E_exciton,"Real","no_rydbergs")
    plot_eigen(imag.(eigenvalues), gap_hBN, E_exciton,"Imag","no_rydbergs")
    print("--> Spectrum plotted! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)

    println(" · Building the laser field...")
    Ttotal, ti, tf = time_prop_param(τ, laser_type)
    time_array = collect(ti:dt:tf)
    println("   · The number of timesteps is ",length(time_array))
    println("   · The total duration of the laser field is ", au_to_fs(Ttotal)," fs")
    E_t = electric_field(time_array, ω, τ, E0, laser_type)
    plot_electric_field(time_array, E_t)
    print("--> Electric field built! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)

    println(" · Constructing the observables...")
    j_operator = current_op(H0,r_operator)
    ground_state = eigenvectors[:,1] * eigenvectors[:,1]' 
    print("--> Observables built! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)

    if projected_bool
        println(" · Starting to construct the system with rydberg states...")
        if optim_bool_ryd == true
            H0_rydbergs, a_rydbergs, η_rydbergs = hamiltonian_optimization(x_array, -gap_hBN, -E_exciton, a_rydbergs_prev, η_rydbergs_prev)
        else
            println(" · We are not optimizing!!!")
            H0_rydbergs = hamiltonian_no_optimization(x_array, a_rydbergs_prev, η_rydbergs_prev)
            a_rydbergs = a_rydbergs_prev
            η_rydbergs = η_rydbergs_prev
        end 
        r_rydbergs = deepcopy(r_operator)
        j_rydbergs = current_op(H0_rydbergs, r_rydbergs)
        print("--> Projection done! "); printstyled(" ✓\n", color = :light_green)
        flush(stdout)
        
        plot_potential(x_array, gap_hBN, E_exciton, a_rydbergs, η_rydbergs, a_no_rydbergs, η_no_rydbergs)

        println(" · Getting the spectrum of the rydberg hamiltonian...")
        #eigenvalues_ryd, eigenvectors_ryd = eigs(H0_rydbergs, nev = 1000, which = :SR)
        eigenvalues_ryd, eigenvectors_ryd = eigen(Array(H0_rydbergs))
        plot_eigen(real.(eigenvalues_ryd), gap_hBN, E_exciton, "Real", "rydbergs")
        plot_eigen(imag.(eigenvalues_ryd), gap_hBN, E_exciton, "Imag", "rydbergs")
        print("--> Spectrum of the rydberg hamiltonian obtained! "); printstyled(" ✓\n", color = :light_green)
        ground_state_rydbergs = eigenvectors_ryd[:,1] * eigenvectors_ryd[:,1]'
        flush(stdout)
    end 

    if php_calc_bool
        println(" · This calculation is also doing the PHP procedure!!")
        println(" · Obtaining the spectrum...")
        H_php, a_php = hamiltonian_optimization_ground(x_array, -gap_hBN, 1.4142)
        eigenvalues_php, eigenvectors_php = eigen(Array(H_php))
        plot_eigen(real.(eigenvalues_php), gap_hBN, E_exciton, "Real", "php_pre_projection")
        plot_eigen(imag.(eigenvalues_php), gap_hBN, E_exciton, "Imag", "php_pre_projection")
        print("--> Spectrum of the php hamiltonian obtained! "); printstyled(" ✓\n", color = :light_green)
        println(" · Building the P operator")
        j_php = current_op(H_php, r_operator)
        ground_state_php = eigenvectors_php[:,1] * eigenvectors_php[:,1]'
        P = rydbergs_projector(eigenvalues_php, eigenvectors_php)
        println("--> Done!!")
    end
    print(" Static part of the calculation is now complete! "); printstyled(" ✓ ✓\n", color = :light_green)
    flush(stdout)

    if timeprop_bool
        println("-"^50)
        println("Starting the time propagation!!")
        t1 = time()
        initial_state  = eigenvectors[:,1]
        laser_func = interpolate_laser_diffeq(time_array, E_t)
        results = sesolve(H0, [r_operator], [laser_func], initial_state, time_array, [j_operator,r_operator, I, ground_state])
        println("-"^50)
        println("Ending the time propagation!!")
        println("It took a total of $(time()-t1) seconds")
        flush(stdout)
        
        println("Saving the results and plotting them...")
        plot_results(time_array, results, [ω],"no_rydbergs")
        print("--> Done! "); printstyled(" ✓\n", color = :light_green)

        if projected_bool
            println("-"^50)
            println("Starting the time propagation with the presence of rydberg states!!")
            t2 = time()
            results_rydbergs = sesolve(H0_rydbergs, [r_rydbergs], [laser_func], eigenvectors_ryd[:,1], time_array, [j_rydbergs, r_rydbergs, I, ground_state_rydbergs])
            println("-"^50)
            println("Ending the time propagation!!")
            println("It took a total of $(time()-t2) seconds")
            flush(stdout)

            println("Saving the results and plotting them...")
            plot_results(time_array, results_rydbergs, [ω], "rydbergs")
            print("--> Calculation Done! "); printstyled(" ✓\n", color = :light_green)
            flush(stdout)
        end

        if php_calc_bool
            println("-"^50)
            println("Starting the time propagation within the PHP procedure")
            t2 = time()
            results_php = sesolve_projected(H_php, [r_operator], [laser_func], eigenvectors_php[:,1], time_array, [j_php, r_operator, I, ground_state_php],P)
            println("-"^50)
            println("Ending the time propagation!!")
            println("It took a total of $(time()-t2) seconds")
            flush(stdout)

            println("Saving the results and plotting them...")
            plot_results(time_array, results_php, [ω], "php")
            print("--> Calculation Done! "); printstyled(" ✓\n", color = :light_green)
            flush(stdout)
        end
    else
        print("--> Calculation Done! "); printstyled(" ✓\n", color = :light_green)
        flush(stdout)
    end    
end