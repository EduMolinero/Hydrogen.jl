
"""
    Vatom()
Returns a soft atomic potential for the 1d atom, it depends on two parameters, namely, a and η.
a will modulate the depth of the potential while η will modulate the width. 
Playing with these two parameters will help us fix the ground state energy and the energy of the first excited state.
"""
@inline function Vatom(x,a,η)
    return -η/sqrt(x^2 + a^2)
end


"""
    hamiltonian_time_in()
It creates the time-indepedent part of the hamiltonian as the sum of the kinetic energy 
and the potential energy.
The kinetic energy matrix is defined as (for ix > 1):
        T[ix,ix-1] = -0.5 * (1.0/dx^2)
        T[ix-1,ix] = -0.5 * (1.0/dx^2)
        T[ix,ix] = 1.0/dx^2
and the potential matrix as :
    Vatom[ix,ix] = Vatom(xs[ix], a, η)
"""
function hamiltonian_time_in(xs::AbstractArray, a::S, η::S) where S <: Number
    dim = length(xs)
    dx = xs[2] - xs[1]
    t_hopping = -0.5 * (1.0/dx^2)
    t_onsite = 1.0/dx^2 
    ## we start with the diagonal terms
    Isparse_diag = collect(1:dim)
    Jsparse_diag = collect(1:dim)
    Ksparse_diag = (Vatom.(xs, a, η) .+ t_onsite) .* ones(ComplexF64, length(Isparse_diag))

    ### Now we do the tridigonal terms (up and down)
    # We start at index 2 becuase the term H[1,1] has no tridigonal terms
    Isparse_up = collect(2:dim)
    Jsparse_up = Isparse_up .- 1
    Ksparse_up = t_hopping .* ones(ComplexF64,length(Isparse_up))

    Jsparse_down = collect(2:dim)
    Isparse_down = Jsparse_down .- 1
    Ksparse_down = t_hopping .* ones(ComplexF64,length(Jsparse_down))

    #Append everything together
    append!(Isparse_diag,Isparse_up)
    append!(Isparse_diag,Isparse_down)
    
    append!(Jsparse_diag,Jsparse_up)
    append!(Jsparse_diag,Jsparse_down)

    append!(Ksparse_diag,Ksparse_up)
    append!(Ksparse_diag,Ksparse_down)    

    return sparse(Isparse_diag, Jsparse_diag, Ksparse_diag, dim,dim)
end

"""
    min_energy()
Returns the following quantity:
   L[Δ, E_g, a, η] = (E_0[a, η] - Δ)^2 + (E_1[a, η] - E_g)^2
This will help us to fix the energy of the ground state and the excited state.
"""
function min_energy(param, ground_state_desired, first_state_desired, x)
    H0 = hamiltonian_time_in(x,param[1],param[2])
    eigvals, eigvecs = eigs(H0, nev = 6, which = :SR, maxiter = 20000, tol = 1e-4)
    return abs2(ground_state_desired - eigvals[1]) + abs2(first_state_desired - eigvals[2])
end


"""
    fix_states()
This functions return us
        min_{a, η} L[Δ, E_g, a, η]
,i.e., it gives for which values of (a, η) the ground state energy is equal to Δ and also the first excited state is equal to E_exciton.
This will help us model the excitons in 2d crystals as states of 1d atoms. 
a0, η0 are the starting parameters which commonly will be a0=1.4142 and η0 = 1.0
"""
function fix_states(x::AbstractArray, Δ::S, E_exciton::S, a0::S, η0::S) where S <: Number
    res = optimize(param -> min_energy(param, Δ, E_exciton, x), [a0, η0],
                Optim.Options(x_tol=1e-2, f_tol = 1e-2, show_trace=true))
    flush(stdout)
    return Optim.minimizer(res)[1], Optim.minimizer(res)[2]
end

"""
    hamiltonian_optimization()
Returns the 1d hamiltonian with ground state energy of Δ and first excited state energy of E_exciton.
This function does the optimization procedure.
"""
function hamiltonian_optimization(x::AbstractArray, Δ::S, E_exciton::S, a0::S, η0::S) where S <: Number 
    println("    ·  Starting minimization procedure...")
    println("    ·  The parameters to choose are ground state energy: ",Δ, " and the first state energy is ", E_exciton)
    flush(stdout)
    a_fixed, η_fixed = fix_states(x, Δ, E_exciton, a0, η0)
    print("    ·  Done! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)
    return hamiltonian_time_in(x, a_fixed, η_fixed), a_fixed, η_fixed
end

"""
    hamiltonian_no_optimization()
Returns the 1d hamiltonian with ground state energy of Δ and first excited state energy of E_exciton.
This function does NOT do the optimization procedure but it depends on the given parameters a and η from previous calculations. 
"""
function hamiltonian_no_optimization(x::AbstractArray, a::S, η::S) where S <: Number 
    println("    ·  Building hamiltonian...")
    return hamiltonian_time_in(x, a, η)
    print("    ·  Done! "); printstyled(" ✓\n", color = :light_green)
end

function R_operator(x::AbstractArray)
    """ Position operator """
    dim = length(x)
    Isparse = collect(1:dim)
    Jsparse = collect(1:dim)
    Vsparse = x .* ones(ComplexF64, dim)
    return sparse(Isparse,Jsparse,Vsparse,dim,dim)
end

function commutator(A::T, B::T) where T <: AbstractArray
    """ Returns   (A*B - B*A)"""
    result = zeros(eltype(A), size(A))
    mul!(result, A, B, 1.0, 0.0 )
    mul!(result, B, A, -1.0, 1.0)
    return result
end

@inline function current_op(H0::T, r_operator::T) where T <: AbstractArray
    """
    Compute the conmumator of the time-indepedent part j = -ev = i * [r,H]
    Note: the time-dependent part does not contribute because H_t ~ r·E ---> [r, r·E] = 0
    """
    return 1im .* commutator(r_operator, H0)
end


function projector_eigenstates!(eigenvalues, eigenstates, projector)
   Threads.@threads for index in eachindex(eigenvalues)
       projector .+=  eigenstates[:, index] * eigenstates[:, index]'  
   end
   return nothing
end

function get_rid_rydbergs(eigenvalues::AbstractArray, eigenvectors::AbstractArray)
   eigenvectors_new = deepcopy(eigenvectors)
   ground_state_energy = minimum(eigenvalues)
   println("   Getting rid of rydbergs states...")
   N_rydbergs = 0
   for (index, eig) in enumerate(eigenvalues)    
       if (ground_state_energy < eig < 0.0)
           N_rydbergs += 1
           eigenvectors_new[:,index] .*= (0.0 + 0.0im) 
       end
   end 
   println("   A total number of ", N_rydbergs, " Rydbergs states have been setted to zero!")
   return eigenvectors_new
end

"""
    rydbergs_projector()
In order to model a system with no excitons (i.e. no interactions), 
we will get rid of the bound states (except the ground states). This is achieved by constructing the projector P 
to the subspace with no bound states.
We build the projector as P = I - Σ_{rydbergs} |n><n| 
"""
function rydbergs_projector(eigenvalues::AbstractArray, eigenvectors::AbstractArray)
   P = zeros(ComplexF64, (length(eigenvalues), length(eigenvalues)))
   eigenvectors_new = get_rid_rydbergs(eigenvalues, eigenvectors) 
   projector_eigenstates!(eigenvalues, eigenvectors_new, P)
   #P = I - Q
   return (P+P')./2
end 

"""
    min_energy_ground()
Returns the following quantity:
   L[Δ, a] = (E_0[a] - Δ)^2 
This will help us to fix the energy of the ground state.
"""
function min_energy_ground(param, ground_state_desired, x)
    H0 = hamiltonian_time_in(x,param[1],1.0)
    eigvals, eigvecs = eigs(H0, nev = 6, which = :SR, maxiter = 20000, tol = 1e-4)
    return abs2(ground_state_desired - eigvals[1])
end


"""
    fix_states_ground()
This functions return us
        min_{a} L[Δ, a]
,i.e., it gives for which values of a the ground state energy is equal to Δ .
This will help us model the excitons in 2d crystals as states of 1d atoms. 
a0, η0 are the starting parameters which commonly will be a0=1.4142.
"""
function fix_states_ground(x::AbstractArray, Δ::S, a0::S) where S <: Number
    res = optimize(param -> min_energy_ground(param, Δ, x), [a0],
                Optim.Options(x_tol=1e-2, f_tol = 1e-2, show_trace=true))
    flush(stdout)
    return Optim.minimizer(res)[1]
end

"""
    hamiltonian_optimization_ground()
Returns the 1d hamiltonian with ground state energy of Δ.
This function does the optimization procedure.
"""
function hamiltonian_optimization_ground(x::AbstractArray, Δ::S, a0::S) where S <: Number 
    println("    ·  Starting minimization procedure...")
    println("    ·  The parameters to choose are ground state energy: ", Δ)
    flush(stdout)
    a_fixed = fix_states_ground(x, Δ, a0)
    print("    ·  Done! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)
    return hamiltonian_time_in(x, a_fixed, 1.0), a_fixed 
end




