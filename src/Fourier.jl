function fourier_transform(array_time::T, array_f::T) where T <: AbstractArray
    """
    Computes the FFT of the given array array_f and time-domain array_time.
    Array in time must be t[i]-t[i-1] = constant
    The FT is defined as: FT = int (exp(-i*w*t) f(t) dt)
    """
    #Some zero padding at the end to avoid noise
    N=length(array_time)
    NFFT=1
    i=0
    while N > NFFT
        i += 1
        NFFT=2^i
    end
    dt=array_time[2] - array_time[1]
    if (NFFT!=N)
        array_time = append!(array_time, array_time[end] .+ collect(range(1,NFFT - N,step=1)).* dt .+ dt )
        array_f = append!(array_f, zeros(ComplexF64, NFFT-N))
    end

    # frequency domain
    t0=array_time[1]
    total_time=dt * N
    dω = (2.0*π)/total_time
    array_ω = LinRange(0, N*dω, N) |> collect

    #FFT
    array_FT_f = fft(array_f)
    array_FT_f .*= dt
    array_FT_f = array_FT_f[1:length(array_ω)] 

    #### Just the fact that the array starts at t0
    if t0 != 0.0
        shift_factor = exp.(-1.0im * t0 .* array_ω)
        array_FT_f .*= shift_factor
    end
    return array_ω, array_FT_f
end 

function envelope_sin2(ti::T, tf::T ,Σ::T, t::T) where T <: Real
    """
    Envolope function for the time-window fourier transform, see Eq.5.2.45 of Rui's thesis.
    This is done to avoid numerical artifacts. 
    """
    if t <= ti
        return 0.0
    end
    if t >= tf
        return 0.0
    end
    if t>=ti+Σ && t<=tf-Σ
        return 1.0
    end
    if t>ti && t<ti+Σ
        temp = sin((π * (t-ti)) / (2.0 * Σ))
        return temp*temp
    end
    if t>tf-Σ && t<tf
        temp = sin((π * (tf-t)) / (2.0 * Σ))
        return temp*temp
    end
end

function envelope_hhg(ti::T, tf::T, Σ::T, array_t::AbstractArray) where T <: Real
    array_env = zeros(eltype(array_t), length(array_t))
    for (index,t) in enumerate(array_t)
        array_env[index] = envelope_sin2(ti,tf,Σ,t)
    end
    return array_env
end

function harmonic_analysis(time_raw::AbstractArray, array_f_raw::AbstractArray, Omega::Real, OmegaFactor = 2)
    """
    This function computes the harmonic analysis of the current j(t).
    It plots the harmonic spectrum up to cutoff frequency.
    """
    
    #FFT only works with evenly spaced arrays, so we interpolate the raw data to fulfill that.
    interpol = LinearInterpolation(time_raw, array_f_raw)
    time = LinRange(time_raw[1], time_raw[end], length(time_raw)) |> collect
    array_f = interpol(time)
    array_time = deepcopy(time)

    #setup stuff
    SigmaEnvelope = fs_to_au(4.0)
    MinValue = 1.0e-40
    OmegaMax = Omega * 100
    light_vel_au = 137.036
    prefactor = 2.0/(3.0 * π *(light_vel_au^3))
    time_i = array_time[1]
    time_f = array_time[end]
    envelope_hhg_temp = envelope_hhg(time_i, time_f, SigmaEnvelope, array_time)
    array_f .*= envelope_hhg_temp

    ## FFT transform
    array_ω, FFT = fourier_transform(array_time,array_f)
    FFT_phase = angle.(FFT)
    FFT = abs.(FFT).^2
    ##### Rescale properly to get proper Harmonic spectrum
    FFT .*= (array_ω .^ OmegaFactor)
    FFT .*= prefactor
    FFT .+= MinValue
    return array_ω[array_ω .< 4.0*OmegaMax], real.(FFT[array_ω .< 4.0*OmegaMax])
end