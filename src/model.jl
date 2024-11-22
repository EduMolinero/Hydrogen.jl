
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
    println("   --> Starting minimization procedure...")
    println("   --> The parameters to choose are ground state energy: ",Δ, " and the first state energy is ", E_exciton)
    flush(stdout)
    a_fixed, η_fixed = fix_states(x, Δ, E_exciton, a0, η0)
    print("   --> Done! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)
    return hamiltonian_time_in(x, a_fixed, η_fixed), a_fixed, η_fixed
end

"""
    hamiltonian_no_optimization()
Returns the 1d hamiltonian with ground state energy of Δ and first excited state energy of E_exciton.
This function does NOT do the optimization procedure but it depends on the given parameters a and η from previous calculations. 
"""
function hamiltonian_no_optimization(x::AbstractArray, a::S, η::S) where S <: Number 
    println("   --> Building hamiltonian...")
    return hamiltonian_time_in(x, a, η)
    print("   --> Done! "); printstyled(" ✓\n", color = :light_green)
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
    println("   --> Starting minimization procedure...")
    println("   --> The parameters to choose are ground state energy: ", Δ)
    flush(stdout)
    a_fixed = fix_states_ground(x, Δ, a0)
    print("   --> Done! "); printstyled(" ✓\n", color = :light_green)
    flush(stdout)
    return hamiltonian_time_in(x, a_fixed, 1.0), a_fixed 
end




###################################################################################################
#############################  OLD VERSIONS OF PROJECTIONS OF RYDBERGS ############################
###################################################################################################



# """
#     rydbergs_index()
# Returns the index of the first (I) and the last (J) rydberg states.
# We will use those indexes to get rid of those states.
# """
# function rydbergs_index(eigenvalues::AbstractArray)
#     ground_state_energy = minimum(eigenvalues)
#     I = findfirst(x -> x > ground_state_energy, eigenvalues)
#     J = findfirst(x -> x > 0.0, eigenvalues)
#     println("   -> This system has ", J-1, " Rydbergs states!!")
#     return I, J
# end

# """
#     diag_rotation_matrix()
# Return the matrix U such that it diagonalize H, i.e. :
#     H_{diag} = U^{dagger} * H * U
# """
# function diag_rotation_matrix(eigenvectors::AbstractArray)
#     dim = length(eigenvectors[:,1])
#     U = eigenvectors[:,1]
#     for index in 2:dim
#         U = hcat(U, eigenvectors[:,index])
#     end
#     return U
# end

# """
#     rotate_to_diag_basis()
# Computes the observable O in the diagonal basis of H.
# """
# function rotate_to_diag_basis(O::AbstractArray, U::AbstractArray)
#     result = zeros(eltype(U), size(O))
#     temp = zeros(eltype(U), size(O))
#     mul!(temp, O, U, 1.0, 0.0)
#     mul!(result, U', temp, 1.0, 0.0)
#     return result
# end

# """
#     reduce_operator()
# Eliminates the elements of the given operator O in the diagonal basis between indexed I and J of the spectrum.
# This function help us projected out the states related to the rydberg states. 
# Therefore, they won't be present during the time propagation.
# """
# function reduce_operator(O::AbstractArray, U::AbstractArray, I::T, J::T) where T <: Number
#     O_diag = rotate_to_diag_basis(O, U)
#     reduced_O = O_diag[[collect(1:I-1); collect(J+1:end)], [collect(1:I-1); collect(J+1:end)]]
#     return reduced_O 
# end    

# """
#     reduce_U()
# This function reduce the matrix for the basis change U as in the case of reduce_operator().
# The purpose of this is to be able to be back to original basis after we eliminate the rydberg states.
# """
# @inline function reduce_U(U::AbstractArray, I::T, J::T) where T <: Number
#     return U[[collect(1:I-1); collect(J+1:end)], [collect(1:I-1); collect(J+1:end)]]
# end 

# """
#     project_O_rydberg()
# Main function to the projection scheme. It computes such projection over a list of observables.
# It does the following:
#     1) Get the index of the rydberg states in the diagonal basis of H
#     2) Build the matrices that do the change of basis to the original basis of H to its diagonal basis.
#     3) Rotate all observables to such basis and then, eliminate the elements related to the rydbergs states.
#     4) Rotate the reduced O back to the original basis to perform the time propagation scheme.
# """
# function project_O_rydberg(O::AbstractArray, eigenvalues::AbstractArray, eigenvectors::AbstractArray)
#     # Get the rydberg index
#     I, J = rydbergs_index(eigenvalues)
#     # Rotation matrix
#     println("   -> Building U matrices....")
#     U = diag_rotation_matrix(eigenvectors)
#     #U_reduced = reduce_U(U, I, J)
#     print("   -> Done!!!"); printstyled(" ✓\n", color = :light_green)
#     println("   -> Projecting one observable....")
#     O_reduced_diagonal = reduce_operator(O, U, I, J)
#     print("   -> Done!!!"); printstyled(" ✓\n", color = :light_green)
#     # for O in O_list
#     #     O = Array(O)
#     #     O_reduced_diagonal = reduce_operator(O, U, I, J)
#     #     #O_reduced_nondiagonal = rotate_to_diag_basis(O_reduced_diagonal, U_reduced)
#     #     hcat(O_list_reduced, O_reduced_diagonal)
#     #     @show size(O)
#     #     @show size(O_reduced_diagonal)
#     #     #@show size(O_reduced_nondiagonal)

#     # end
#     # @show typeof(O_list_reduced)
#     # @show size(O_list_reduced)
#     return O_reduced_diagonal
# end

#function projector_eigenstates!(eigenvalues, eigenstates, projector)
#    Threads.@threads for index in eachindex(eigenvalues)
#        projector .+=  eigenstates[:, index] * eigenstates[:, index]'  
#    end
#    return nothing
#end
#
#function get_rid_rydbergs(eigenvalues::AbstractArray, eigenvectors::AbstractArray)
#    eigenvectors_new = deepcopy(eigenvectors)
#    ground_state_energy = minimum(eigenvalues)
#    println("   Getting rid of rydbergs states...")
#    N_rydbergs = 0
#    for (index, eig) in enumerate(eigenvalues)    
#        if (ground_state_energy < eig < 0.0)
#            N_rydbergs += 1
#            eigenvectors_new[:,index] .*= (0.0 + 0.0im) 
#        end
#    end 
#    println("   A total number of ",N_rydbergs, " Rydbergs states have been setted to zero!")
#    return eigenvectors_new
#end
#
#function rydbergs_projector(eigenvalues::AbstractArray, eigenvectors::AbstractArray)
#    """
#    In order to model a system with no excitons (i.e. no interactions), 
#    we will get rid of the bound states (except the ground states). This is achieved by constructing the projector P 
#    to the subspace with no bound states.
#    We build the projector as P = I - Σ_{rydbergs} |n><n| 
#    """
#    P = zeros(ComplexF64, (length(eigenvalues), length(eigenvalues)))
#    eigenvectors_new = get_rid_rydbergs(eigenvalues, eigenvectors) 
#    projector_eigenstates!(eigenvalues, eigenvectors_new, P)
#    return (P+P')./2
#end
#
#function build_projector(P::AbstractArray, O::AbstractArray)
#    """ Construct the operator P * O * P from operator O and projector P"""
#    temp_1 = zeros(eltype(O), size(O)) 
#    temp_2 = zeros(eltype(O), size(O)) 
#    mul!(temp_1, O, P, 1.0,  0.0)
#    mul!(temp_2, P, temp_1, 1.0, 0.0)
#    return (temp_2 + temp_2')./2
#end
#