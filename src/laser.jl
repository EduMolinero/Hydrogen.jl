"""
    electric_field()
Returns the electric_field for a laser pulse with formula
        E(t) = E_0 * sin(ωt) * f(t,τ,ω)
where f(t) is some envelope function. Since it is meant for the one-dimensional case, 
we do not need to worry about ellipticity. We assume that CEP is zero and the field is also centered around zero.
"""
function electric_field(time::AbstractArray, ω::T, τ::T, E0::T, type_of_wave::AbstractString) where T <: Real
    η = ω .* time
    envelope = envelope_laser.(η, ω, τ, type_of_wave)
    return E0 .* envelope .* sin.(η)
end

"""
    envelope_laser()
Returns the envelope f(t) for the laser field. Currently, it only supports cos2 and gaussian envelopes.
For more info on the parameters see the book "Atoms in intense laser fields" by C. J. Joachain.
"""
function envelope_laser(η::S, ω::S, τ::S, type_of_wave::AbstractString) where S <: Number
    if type_of_wave == "cos2"
        if abs(η) > π * ω * τ * 0.5 
            return 0.0
        end
        result = cos(η / (ω*τ))^2 
    elseif type_of_wave == "gaussian"
        temp = - 1.0 * (η^2) / (2.0 * (ω * τ)^2)
        result =  exp(temp)
    end
    return result
end 

"""
    time_prop_param()
Returns some laser parameters depending on which envelope we have decided to use.
"""
function time_prop_param(τ::Real, type_of_wave::AbstractString)
    if type_of_wave == "cos2"
        Ttotal = π * τ
    elseif type_of_wave == "gaussian"
        Ttotal = 8.0 * τ
    end
    ti = - Ttotal * 0.5
    tf = Ttotal * 0.5
    return Ttotal, ti, tf 
end


@inline function interpolate_laser_diffeq(time::AbstractArray, laser::AbstractArray)
    return LinearInterpolation(time, laser)
end