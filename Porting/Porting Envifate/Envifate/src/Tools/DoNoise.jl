module DoNoise
"""
"""



function transmission_loss( r::Float64 )
    return 20 * log10(r)
end



function atmospheric_absorpion_loss( r::Float64, height_m::Float64, relative_humidity::Float64, temperature_k::Float64, frequency::Float64 )
    # Calculate atmospheric absorption coefficient using ANSI S1.26-1995 standard

    # Convert elevation to atmospheric pressure
    atmospheric_pressure = 101.325 * ( 1 - ( 2.25577 * ( 10^(-5) ) * height_m ) )^5.25588

    # Convert relative humidity to molar concentration of water vapor
    C = ( -6.8346 * ( ( 273.16 / temperature_k )^1.261 ) ) + 4.6151
    p_saturation_pressure = 10^C
    humidity = relative_humidity * p_saturation_pressure * ( ( atmospheric_pressure / 101.325 )^(-1) )
    # Calculate derived values for subsequent equations
    p_atm_pressure = atmospheric_pressure / 101.325
    T_Tr = temperature_k / 293.15

    # Calculate frO (equation 3)
    O_relax_freq = ( p_atm_pressure * ( ( 24 + ( 4.04 * (10^4) ) ) * humidity ) * (0.02 + humidity) ) / ( 0.391 + humidity )

    # Calculate frN (equation 4)
    N_relax_freq = p_atm_pressure * ( T_Tr^(-0.5) ) * ( 9 + ( 280 * humidity * ( ℯ^( -4.170 * ( ( T_Tr^(-0.33333) ) - 1 ) ) ) ) )

    # Calculate alpha (equation 5)
    term1 = 1.84 * (10^(-11)) * (p_atm_pressure^(-1)) * (T_Tr^0.5)
    term2 = (T_Tr^(-2.5)) * ( 0.01275 * (ℯ^( -2239.1 / temp_k )) * ( O_relax_freq / ( (O_relax_freq^2) + (freq^2) ) ) )
    term3 = 0.1068 * (ℯ^( -3352 / temp_k )) * ( N_relax_freq / ( (N_relax_freq^2) + (freq^2) ) )
    
    α = 8.686 * (frequency^2) * ( term1 + term2 + term3 )

    return α * r / 100
end

function ground_loss( ground::Matrix, source_xy::Tuple{Float64, Float64} )

end

function vegetation_loss( vegetation_map::Matrix, source_xy::Tuple{Float64, Float64} )
end

function wind_loss()
end



end # module