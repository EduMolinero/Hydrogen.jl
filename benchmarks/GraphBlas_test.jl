using Optim
using LinearAlgebra, SparseArrays, Arpack
using BenchmarkTools
using SuiteSparseGraphBLAS  


function op_sum_before_rk(A, B,vec)
    M = A + B
    temp1 = M * vec
    temp2 = M * (vec + 0.5 * temp1)
    temp3 = M * (vec + 0.5 * temp2)
    temp4 = M * (vec + temp3)
    return vec + (1/6) * (temp1 + 2 * temp2 + 2 * temp3 + temp4)
end

function op_sum_after_rk(A, B,vec)
    temp1_a = A * vec
    temp2_a = A * (vec + 0.5 * temp1_a)
    temp3_a = A * (vec + 0.5 * temp2_a)
    temp4_a = A * (vec + temp3_a)

    temp1_b = B * vec
    temp2_b = B * (vec + 0.5 * temp1_b)
    temp3_b = B * (vec + 0.5 * temp2_b)
    temp4_b = B * (vec + temp3_b)

    result_a = vec + (1/6) * (temp1_a + 2 * temp2_a + 2 * temp3_a + temp4_a)
    result_b = vec + (1/6) * (temp1_b + 2 * temp2_b + 2 * temp3_b + temp4_b)

    return result_a + result_b
end

n = 10000
d = 0.05

A = sprand(n,n,d)
B = sprand(n,n,d)
vec = rand(n)

println("Benchmarking...")
println("n = $n, d = $d")
println("Julia threads: $(Threads.nthreads())")

println("op_sum_before_rk no GraphBLAS")
@btime op_sum_before_rk(A,B,vec)

println("op_sum_after_rk no GraphBLAS")
@btime op_sum_after_rk(A,B,vec)

println("op_sum_before_rk using GraphBLAS")
@btime op_sum_before_rk(GBMatrix(A), GBMatrix(B), GBVector(vec))

println("op_sum_after_rk using GraphBLAS")
@btime op_sum_after_rk(GBMatrix(A), GBMatrix(B), GBVector(vec))

