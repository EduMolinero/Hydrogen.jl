using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

function sesolve(H0,H0_td,fs_H0_td,psi0,ts,e_ops ;dt=nothing,rtol=1e-10,atol=1e-10)
    """
    H0 is the time independent Hamiltonian
    H0_td is a list of a operators such as 
    H(t) = H0 + sum_i H0_td[i] fs_H0_td[i](t)
    psi0 is the initial wavefunction
    ts is an array of times for the propagation
    e_ops is a list of operators to perform expectation values
    """ 

    temp = psi0 .* 0.0
    results = zeros(ComplexF64,length(ts),length(e_ops))
    ps = [H0,H0_td,fs_H0_td,temp]
    if dt === nothing
        dt=ts[2]-ts[1]
    end
    
    function minus_im_Heff_sesolve!(psires,psi,ps,t)
        H0,H0_td,fs_H0_td,temp = ps
        # Hamiltonian part
        mul!(psires,H0,psi,-1.0im,0.0)
        for i in eachindex(H0_td)
            mul!(psires,H0_td[i],psi,-1.0im*fs_H0_td[i](t),1.0)
        end
        return psires
    end 
    
    function get_callback(ts,results,e_ops,temp)
        function savefunc(u,t,integrator)
            it=findall(ts.≈t)[1]
            for i_eop in eachindex(e_ops)
                results[it,i_eop] = expectation_value(u,e_ops[i_eop],temp)
            end
            return nothing
        end
        FunctionCallingCallback(savefunc; funcat=ts)
    end

    cb2=get_callback(ts,results,e_ops,temp)
    cb=CallbackSet([cb2]...);
    prob = ODEProblem(minus_im_Heff_sesolve!,psi0,(ts[1],ts[end]),ps)
    sol = solve(prob,DP5(),save_start=false,save_end=false,save_everystep=false,callback=cb,
    reltol=rtol,abstol=atol,maxiters=1e6,dt=dt, progress = true, progress_steps = 1)
    
    tempwf = nothing
    return results
end


function sesolve_projected(H0,H0_td,fs_H0_td,psi0,ts,e_ops,projector;dt=nothing,rtol=1e-10,atol=1e-10)
    """
    H0 is the time independent Hamiltonian
    H0_td is a list of a operators such as 
    H(t) = H0 + sum_i H0_td[i] fs_H0_td[i](t)
    psi0 is the initial wavefunction
    ts is an array of times for the propagation
    e_ops is a list of operators to perform expectation values
    """ 

    temp = psi0 .* 0.0
    temp_2 = psi0 .* 0.0
    results = zeros(ComplexF64,length(ts),length(e_ops))
    ps = [H0,H0_td,fs_H0_td,temp]
    if dt === nothing
        dt=ts[2]-ts[1]
    end
    
    function minus_im_Heff_sesolve!(psires,psi,ps,t)
        H0,H0_td,fs_H0_td,temp = ps
        # Hamiltonian part
        mul!(psires,H0,psi,-1.0im,0.0)
        for i in eachindex(H0_td)
            mul!(psires,H0_td[i],psi,-1.0im*fs_H0_td[i](t),1.0)
        end
        return psires
    end 
    
    function get_callback_measurement(ts,results,e_ops,temp)
        function savefunc(u,t,integrator)
            it=findall(ts.≈t)[1]
            for i_eop in eachindex(e_ops)
                results[it,i_eop] = expectation_value(u,e_ops[i_eop],temp)
            end
            return nothing
        end
        FunctionCallingCallback(savefunc; funcat=ts)
    end

    function get_callback_projection(ts, P, temp)
        function savefunc(u,t,integrator)
            mul!(temp,P,u) 
            return nothing
        end
        FunctionCallingCallback(savefunc; funcat=ts)
    end

    cb_measurement = get_callback_measurement(ts,results,e_ops,temp)
    cb_projection = get_callback_projection(ts, projector, temp_2)
    cb=CallbackSet([cb_projection, cb_measurement]...);
    prob = ODEProblem(minus_im_Heff_sesolve!,psi0,(ts[1],ts[end]),ps)
    sol = solve(prob,DP5(),save_start=false,save_end=false,save_everystep=false,callback=cb,
    reltol=rtol*0.1,abstol=atol*0.1,maxiters=1e6, dt=dt, progress = true, progress_steps = 1)
    
    tempwf = nothing
    return results
end

function expectation_value(wf,oper,tempwf)
    """ Calculates <psi|oper|psi> """
    mul!(tempwf,oper,wf)
    return dot(wf,tempwf)
end

