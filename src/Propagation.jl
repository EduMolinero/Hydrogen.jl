module propagation

using ..RungeKutta


    function propagate_timestep(H0::AbstractArray, H_td::AbstractArray, f_td, ψ::AbstractArray, ti::Real, dt::Real, tableau::RungeKutta.RK_Butcher_tableau)
        """
        Propagate the wavefunction a single timestep dt using the Runge-Kutta method
        ψ(t+dt) = ψ(t) + dt * RK(t,ψ(t)) where RK() is the choosen Runge-Kutta method
        Currently implemented: Euler, RK4, DP5, Tsit5
        https://www.johndcook.com/blog/2020/02/13/runge-kutta-methods/

        Since both H0 and H_td are Sparse matrices, we will split the computation of the RK(t,ψ(t))  in two parts:
            temp_H0 ~ H0 * ψ
            temp_H0 ~ H_td * ψ
        and then sum them up to get the final result
            ψ(t+dt) = ψ(t) + dt * (temp1 + temp2)

        This is faster than summing the two matrices and then multiplying by ψ, becuase that process converts them to dense matrices and then back to sparse matrices.
        """


        temp_H0 = 0.0 .* ψ # This ensures is the same type of array: MKL or CUDA
        temp_Ht = 0.0 .* ψ

        ks_vector_H0 = [0.0 .* ψ for i in 1:tableau.s]
        ks_vector_Ht = [0.0 .* ψ for i in 1:tableau.s]

        # Loop over order of the RK method
        for i in 1:RK_tableau.s
            time_i = ti + RK_tableau.c[i]*dt
            f_td_time_i = f_td(time_i) 

            #temp ~ y_n + h * sum_j=1^i-1( a_ij * k_j) 
            temp_H0 .= ψ
            temp_Ht .= ψ

            for j in 1:i-1
                if RK_tableau.A[i,j] != 0.0
                   axpy!(dt * RK_tableau.A[i,j], ks_vector_RK_H0[j], temp_H0)
                   axpy!(dt * RK_tableau.A[i,j], ks_vector_RK_Ht[j], temp_Ht)
                end
            end

            ks_vector_RK_H0[i] .= 0.0
            ks_vector_RK_H0[i] .= 0.0
             
            ## Update temp1 
            # k_i ~ f(t_i, temp )
            mul!(ks_vector_RK[i], -1.0im .* H0, temp1)
            mul!(ks_vector_RK[i], H_td, temp2)

        end

    end

end