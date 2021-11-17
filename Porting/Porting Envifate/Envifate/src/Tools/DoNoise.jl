module DoNoise
"""
"""



using ArchGDAL
using GeoArrays
using GeoRegions
using GeoStats
using Plots
using Shapefile

Base.convert(::Type{Int64}, n::Float64) = Int64(round(n))
Base.:-( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] - y[1], x[2] - y[2] )
Base.:-( x::Vector{T}, y::Tuple{T, T} ) where {T <: Number} = length(x) == length(y) ? [ e1 - e2 for (e1, e2) in zip(x, y) ] : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:-( x::Tuple{T, T}, y::Vector{T} ) where {T <: Number} = length(x) == length(y) ? Tuple( e1 - e2 for (e1, e2) in zip(x, y) ) : throw(ArgumentError("`x` and `y` must have the same size"))
Base.:+( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] + y[1], x[2] + y[2] )
Base.:*( x::Tuple{Number, Number}, y::Number ) = ( x[1] * y, x[2] * y )
Base.:*( x::Number, y::Tuple{Number, Number} ) = y * x
Base.:*( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] * y[1], y[1] * y[2] )
Base.:/( x::Tuple{Number, Number}, y::Number ) = ( x[1] / y, x[2] / y )
Base.:/( x::Number, y::Tuple{Number, Number} ) = y / x
Base.:/( x::Tuple{Number, Number}, y::Tuple{Number, Number} ) = ( x[1] / y[1], x[2] / y[2] )
Base.:^( x::Tuple{Number, Number}, y::Number ) = ( x[1]^y, x[2]^y )

toInt( n::Number, method::Function=round )::Int64 =  n isa Int64 ? n : Int64( method(n) )
# Functions to find the coordinates of the point resulting from the rotation of "(xp, yp)" by a angle "θ" around "(xc, yc)" 
rotate_x( xp, yp, xc, yc, θ ) = toInt( (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) + xc )
rotate_y( xp, yp, xc, yc, θ ) = toInt( (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) + yc )
rotate_point( xp, yp, xc, yc, θ ) = ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )

"""
    transmission_loss( r::Real )

Compute the transmission loss of a noise over `r` distance
"""
function transmission_loss( r::Real )
    return 20 * log10(r)
end



function atmospheric_absorpion_loss( r::Real, height_m::Real, relative_humidity::Real, temperature_k::Real, frequency::Real )
    # Calculate atmospheric absorption coefficient using ANSI S1.26-1995 standard

    # Convert elevation to atmospheric pressure
    atmospheric_pressure = 101.325 * ( 1 - ( 2.25577 * 10^(-5) * height_m ) )^5.25588

    # Calculate derived values for subsequent equations
    p_atm_pressure = atmospheric_pressure / 101.325
    t_tr = temperature_k / 293.15

    # Convert relative humidity to molar concentration of water vapor
    C = ( -6.8346 * ( 273.16 / temperature_k )^1.261 ) + 4.6151
    p_saturation_pressure = 10^C
    humidity = relative_humidity * p_saturation_pressure * p_atm_pressure^(-1)

    # Calculate relaxation frequency of O (equation 3)
    #   frO₂ = ( p_atm_pressure * ( (24 + 4.04e04) * humidity ) * (0.02 + humidity) ) / (0.391 + humidity)
    frO₂ = p_atm_pressure * (
                               24 + ( 
                                        4.04e04 * humidity * (
                                                                 (0.02 + humidity) / (0.391 + humidity)
                                                             )
                                    )
                           )

    # Calculate relaxation frequency of N (equation 4)
    frN₂ = p_atm_pressure * √t_tr * (
                                       9 + (
                                               280 * humidity * ℯ^( -4.170 * (
                                                                                 t_tr^(-1/3) - 1
                                                                             )
                                                                  )
                                           )
                                   )

    # Calculate alpha (equation 5)
    term1 = 1.84*10^(-11) * p_atm_pressure^(-1) * √t_tr
    #   term2 = t_tr^(-2.5) * ( 0.01275 * ℯ^(-2239.1 / temp_k) * ( frO₂ / (frO₂^2 + freq^2) ) )
    term2 = t_tr^(-2.5) * (
                              0.01275 * ℯ^(-2239.1 / temp_k) / (
                                                                   frO₂ + (freq^2 / frO₂)
                                                               )
                          )
    #   term3 = 0.1068 * ℯ^(-3352 / temp_k) * ( frN₂ / (frN₂^2 + freq^2) )
    term3 = t_tr^(-2.5) * (
                              0.1068 * ℯ^(-3352 / temp_k) / (
                                                                frN₂ + (freq^2 / frN₂)
                                                            )
                          ) 
    #   α = 8.686 * (frequency^2) * ( term1 + term2 + term3 )
    α = frequency^2 * (term1 + term2 + term3)

    return α * r / 100
end



function minmax( profile::Vector, rel_h_src::Real=0.0, rel_h_rec::Real=0.0 )::Vector{Tuple{Int64, Float64}}
    if length(profile) <= 1
        return [ (1, -rel_h_src), (1, -rel_h_rec) ]
    end

    # Test for the profile having a finite length; will be zero if directly overhead
    if abs( profile[1][1] - profile[end][1] ) < 10.0
        return [ (1, rel_h_src), (length(profile), rel_h_rec) ]
    end

    # Define parameters A and B such that height of the source-receiver line
    # is given by z = Ax + B.
    x1 = profile[1][1]
    z1 = profile[1][2] + rel_h_src
    xn = profile[end][1]
    zn = profile[end][2] + rel_h_rec
    a = (zn - z1) / (xn - x1)
    b = z1 - a*x1

 #=
    # Look for the max and min
    prf = map( row -> ( row[1], row[2] - ( a*row[1] + b ) ), profile )
    sort!( prf, by=(row) -> row[2] )
    
    return [ prf[1], prf[end] ]
 =# 
    min = ( 1, profile[1][2] )  
    max = ( 1, profile[1][2] )
    for i in 1:length(profile)
        hgt = profile[i][2] - a*profile[i][1] + b
        if hgt > max[2]
            max = (i, hgt)
        end
        if hgt < min[2]
            min = (i, hgt)
        end
    end
    return [min, max]
end

function delbaz( freq::Real, flow_res::Real )::Complex
    dumr = 1.0 + 9.08 * (flow_res/freq)^0.75
    dumi = -11.9 * (flow_res/freq)^0.73
    return complex(dumr, dumi)
end

function subw( a::Real, b::Real, c::Real, d::Real, w::Complex )::Complex
    an = (a^2 - b^2  - d)^2 + (2.0 * a * b)^2
    r_i = c * (a^2 + b^2 + d)
    wr = b * r_i / an + w.re
    wi = a * r_i / an + w.im
    return complex( wr, wi )
end

function ww(t::Complex)::Complex
    a = abs(t.re)
    b = abs(t.im)

    if a <= 3.9 && b <= 3.0
        a1 = cos(2.0 * a * b)
        b1 = sin(2.0 * a * b)
        c1 = ℯ^( -2.0 * b * π / 0.8 ) - cos( 2.0 * a * π / 0.8 )
        d1 = sin( 2.0 * a * π / 0.8 )

        pq = 2.0 * ℯ^( -( a^2 + 2.0 * b * π / 0.8 - b^2 ) ) / ( c1^2 + d1^2 )
        p2 = pq * ( a1*c1 - b1*d1 )
        q2 = pq * ( a1*d1 - b1*c1 )

        aa = a^2
        bb = b^2
        ah = 0.8 * b / ( π * (aa + bb) )
        ak = 0.8 * a / ( π * (aa + bb) )

        for i in 1:5
            an = ( bb - aa + i^2 * 0.64 )^2 + 4.0 * aa * bb
            hk = 2.0 * 0.8 / π * ℯ^( -i^2 * 0.64 ) / an
            ah += b * hk * ( bb + aa + 0.64 * i^2 )
            ak += a * hk * ( bb + aa - 0.64 * i^2 )
        end
    
        if b < π/0.8
            ah += p2
            ak -= q2
        end

        if b == π/0.8
            ah += p2 / 2.0
            ak -= q2 / 2.0
        end

        w = complex(ah, ak)
    else
        cd = [
                 (0.4613135, 0.1901635),
                 (0.09999216, 1.7844927),
                 (0.002883894, 5.5253437)
             ]
        w = 0.0 + 0.0im
        for (c, d) in cd
            w = subw(a, b, c, d, w)    
        end
    end

    if t.re < 0.0
        w = conj(w)
        a = -a
    end

    if t.im < 0.0
        wr = 2.0 * ℯ^(b^2 - a^2) * cos(2*a*b) - w.re
        wi = 2.0 * ℯ^(b^2 - a^2) * sin(2*a*b) - w.im
        w = complex(wr, wi)
    end

    return w
end

function qq( r::Real, h::Real, waveno::Real, z::Complex )::Complex
    c = abs(h) / r
    n = (z.re * c + 1.0)^2 + (z.im * c)^2
    
    rr = ( ( c * abs2(z) ) - 1.0 ) / n
    ri = 2.0 * z.im * c / n

    nyr = z.re / abs2(z)
    nyi = -z.im / abs2(z)

    dumr = √(waveno * r / 4.0) * (nyr + c - nyi)
    dumi = √(waveno * r / 4.0) * (nyr + c + nyi)
    t = complex(dumr, dumi)

    w = ww(t)

    fr = 1.0 + √π * -imag(t*w)
    fi = √π * real(t*w)

    dumr = rr + (1.0 - rr) * fr + fi * ri
    dumi = ri + fi * (1.0 - rr) - ri * fr
    return z.re >= 1000.0 ? 1.0+0im : complex(dumr, dumi)
end

function qq2( d::Real, src_h::Real, rec_h::Real, freq::Real, z::Complex )::Complex
    waveno = 2.0 * π * freq / 340.0
    #   r1 = √( d^2 + (src_h - rec_h)^2 )
    r = √( d^2 + (src_h + rec_h)^2 )
    c = (src_h + rec_h) / r
    
    n = (z.re * c + 1.0)^2 + (z.im * c)^2
    rr = ( ( c * abs(z) )^2 - 1.0 ) / n
    ri = 2.0 * z.im * c / n

    nyr = z.re / abs2(z)
    nyi = -z.im / abs2(z)

    dumr = √(waveno * r / 4.0) * (nyr + c - nyi)
    dumi = √(waveno * r / 4.0) * (-nyr - c + nyi)
    t = complex(dumr, dumi)

    w = ww(-t)

    fr = 1.0 + √π * imag(t*w)
    fi = -√π * real(t*w)

    dumr = rr + (1.0 - rr) * fr + fi * ri
    dumi = ri + fi * (1.0 - rr) - ri * fr
    return complex(dumr, dumi)
end

function egal( d1::Real, d2::Real, src_h::Real, rec_h::Real, src_flow_res::Real, rec_flow_res::Real, e_wind_vel::Real, transition_height::Real, turbulence::Real, freq::Real, dum::Real )
 #=
    arg1  = d1
    arg2  = d2
    arg3  = src_h
    arg4  = rec_h
    arg5  = src_flow_res
    arg6  = rec_flow_res
    arg7  = e_wind_vel
    arg8  = transition_height
    arg9  = turbulence
    arg10 = freq
    arg11 = dum
 =#
    src_rx = delbaz( freq, src_flow_res )
    rec_rx = delbaz( freq, rec_flow_res )

    rd = √( (d1 + d2)^2 + (src_h - rec_h)^2 )
    rr = √( (d1 + d2)^2 + (src_h + rec_h)^2 )

    la = 340.0 / freq
    fixed_speed = 2.0 * π * la^(-1)
    e_turbulence_scale = 0.0
    ar = 0.0

    kv = [ 0.0, 0.0, 0.0 ]
    db = [ 0.0+0.0im, 0.0+0.0im, 0.0+0.0im ]

    if e_wind_vel != 0
 #   NON SO SE VADE USATO `turbulence` O COSA
        m = round( transition_height / ( la / 6.0 ) )
        fixed_speed = 2.0 * π * freq / (
                                 340.0 + (
                                             ( m-1 * la/6.0 + la/10.0 ) / 2.0 
                                         ) * e_wind_vel/10.0
                             )
        e_turbulence_scale = 10e3 / freq
    end

    rr = [
             √( d1^2 + (src_h - ha)^2 ),
             √( d1^2 + (src_h + ha)^2 ),
             √( d2^2 + (ha - rec_h)^2 ),
             √( d2^2 + (ha + rec_h)^2 ) 
         ]

    eq(x, rrn, rrm ) = ℯ^complex( 0.0, -x*(rrn + rrm) )
    eq_cr( rrn, rrm ) = rrm * √( rrm * rrn * (rrn + rrm) )
    eq_jarr( kv, k, rrn, rrm; q1=nothing, q2=nothing ) =  ( eq(kv, rrn, rrm) - eq(k, rrn, rrm) ) * ( isnothing(q1) ? 1 : q1 ) * ( isnothing(q2) ? 1 : q2 ) / complex( eq_cr(rrn, rrm), 0.0 ) 

    for j in 1:m
        ha = (j-1) * (la/6.0) + (la/10.0)
 #   NON CAPISCO SE CI VA `turbulence` O COSA
        v = [ e_wind_vel + turbulence / 10.0 * cos( ha * 2.0 * π / e_turbulence_scale ) ]
        push!( 
            v,
            2.0 * e_wind_vel - v[1],
            e_wind_vel
        )

        for i in 1:3
            kv[i] = 2.0 * π * freq / ( 340.0 + (ha/2.0) * v[i]/10.0 )

            #   cr = eq_cr( rr[1], rr[3] )
            #   jarray[i] = ( eq( kv[i], rr[1], rr[3] ) - eq( fixed_speed, rr[1], rr[3] ) ) / (cr + 0.0im)
            #   db[i] += jarray[i]
  
            #   q2 = qq2( d2, ha, rec_h, freq, rec_rx )
            #   cr = eq_cr( rr[1], rr[4] )
            #   jarray[i] = ( eq( kv[i], rr[1], rr[4] ) - eq( fixed_speed, rr[1], rr[4] ) ) * q2 / (cr + 0.0im)
            #   db[i] += jarray[i]
            #   
            #   q1 = qq2( d1, src_h, ha, freq, src_rx )
            #   cr = eq_cr( rr[2], rr[3] )
            #   jarray[i] = ( eq( kv[i], rr[2], rr[3] ) - eq( fixed_speed, rr[2], rr[3] ) ) * q1 / (cr + 0im)
            #   db[i] += jarray[i]
  
            #   cr = eq_cr( rr[2], rr[4] )
            #   jarray[i] = ( eq( kv[i], rr[2], rr[4] ) - eq( fixed_speed, rr[2], rr[4] ) ) * q1 * q2 / (cr + 0im)
            #   db[i] += jarray[i]
            jarray = [
                eq_jarr( kv[i], fixed_speed, rr[1], rr[3] ),
                eq_jarr( kv[i], fixed_speed, rr[1], rr[4], q2=q2 ),
                eq_jarr( kv[i], fixed_speed, rr[2], rr[3], q1=q1 ),
                eq_jarr( kv[i], fixed_speed, rr[2], rr[4], q1=q1, q2=q2 )  
            ]
            for val in jarray
                db[i] += val
            end
        end
    end

    c = ℯ^complex( 0.0, π/4.0 )
    crh = (la / 6.0) * d2 * √(8.0 * π * k) / (16.0 * π^2)
    c *= crh + 0.0im

    db .*= c

    q2 = qq2( d1+d2, src_h, rec_h, freq, rec_rx )
    c = ℯ^complex( 0.0, -k*rr ) / complex( rr, 0.0 ) * q2
    ch = ℯ^complex( 0.0, -k*rd ) / complex( rd, 0.0 )
    c += ch

    c = c / complex( 4.0*π, 0.0 )

    nff = 16.0 * π * 2.0 * rd^2
    if e_wind_vel != 0
        db .+= c
        ar = sum( x -> abs(x)^2, db )
        levturb = 4.3429 * log(ar * nff / 3.0)
        lnot = 4.3429 * log( abs(db[3])^2 * nff )
    else
        levturb = 0.0
        lnot = 4.3429 * log( abs(c)^2 * nff )
    end

    return (levturb, lnot)
end

function varysurf( dists::AbstractVector, ground_type::AbstractVector, src_h::Real, rec_h::Real, soft_atten::Real, hard_atten::Real )::Real
    drefl = src_h * dists[end] / (src_h + rec_h)

    srcInf = ground_type[1] == Soft ?
             src_h > 3.0 ?
                 (20.0 * src_h / 3.0 - 10.0) * src_h :
                 10.0 * src_h :
             30.0 * src_h
       
    recInf = ground_type[end] == Soft ?
             rec_h > 3.0 ?
                 (20.0 * rec_h / 3.0 - 10.0) * rec_h :
                 10.0 * rec_h :
             30.0 * rec_h
    
    srcInf = min( srcInf, 0.7*drefl )
    recInf = min( recInf, 0.7*(dists[end]-drefl) )
    recInf = dists[end] - recInf

    Δl = sum = denom = 0
    for i in 1:length(dists)
        w1 = w2 = 0
        if dists[i] > srcInf && dists[i] < recInf
            w1 = dists[i] < drefl ? rec_h : src_h
            w2, Δl = ground_type[i] == Hard ? (0.5, hard_atten) : (1.0, soft_atten)
        end
        sum += w1 * w2 * Δl
        denom += w1 * w2
    end

    return sum / denom
end

function fres( y::Real )::Complex
    c = 0.797885
    x = c * y
    f = (1.0 + 0.962x) / (2.0 + 1.792x + 3.104x^2 )
    g = 1.0 / (2.0 + 4.142x + 3.492x^2 + 6.67x^3 )
    si = sin( (x / c)^2 )
    co = cos( (x / c)^2 )
    return complex( (-f * si + g * co)/c, (f * co + g * si)/c )
end

function diffraction( r1::Real, a::Real, al2::Real, pm::Real, any::Real, k::Real )::Complex
    df = -ℯ^complex(0.0, k*r1+π/4.0) / complex(r1, 0.0)
    tangent = tan( (π + pm * al2) / (2.0 * any) )

    if tangent != 0
        aa = 1.0 / tangent / (2.0 * any) / √(2.0 * π * k * a)
    else
        # NON SO COSA RAPPRESENTI "aalast" 
        aa = aalast
    end
    aalast = aa
    df *= complex( aa, 0.0 )

    n = pm < -0.9 ?
            al2 > π+any*π ? 1 :
                al2 > π-any*π ? 0 : -1 :
            al2 > any*π-π ? 1 : 0
    xv = 2.0 * k * a * cos( (2.0 * n * any * π - al2) / 2.0 )^2

    aa = -2.0 * √xv
    y = ℯ^complex(0.0, -xv) * fres(√xv)

    return df * y * complex(0.0, aa)
end

function calc_mirror( locs, points; source::Bool )
    i, j, l = source ? (3,2,1) : (3,4,5) 

    # Distance from image to top of hill
    Δx, Δy = source ? points[i]-locs : ( locs[1]-points[i][1], points[i][2]-locs[2] )
    rr = √( Δx^2 + Δy^2 )
    # Angle from image to top of hill
    θi = atan( Δy, Δx )
    # Angle of the hillside
    Δxh, Δyh = source ? points[i]-points[j] : ( points[j][1]-points[i][1], points[i][2]-points[j][2] )
    θh = atan( Δyh, Δxh )

    if θi < θh
        # Angle of the flat
        Δxf, Δyf = source ? points[j]-points[l] : ( points[l][1]-points[j][1], points[j][2]-points[l][2] )
        θf = atan( Δyf, Δxf )
        # Angle of the image path relative to the flat
        θ = θi - θf
        # Net image source-receiver height, as needed by QQ
        Δz = rr * sin(θ)
        rr = abs( rr * cos(θ) )
        return (rr, Δz)
    else
        return nothing
    end
end

function bakkernn( src_loc::AbstractVector, rec_loc::AbstractVector, hills::AbstractVector, src_flow_res::Real, ber_flow_res::Real, rec_flow_res::Real, freq::Real )::Real
    
    # Delany-Bazley under source, berm and reciver
    dbs = delbaz.( Ref(freq), [ src_flow_res, ber_flow_res, rec_flow_res ] )
    waveno = 2.0 * π * freq / 340.0

 # Mirror Source
    #   # Distance from image source to top of hill
     #   Δx, Δy = hills[3] - src_loc[2]
     #   rr = √( Δx^2 + Δy^2 )
     #   # Angle from image to top of hill
     #   θi = atan( Δy, Δx )
     #   # Angle of the hillside
     #   Δxh, Δyh = hills[3] - hills[2]
     #   θh = atan( Δyh, Δxh )
     #   
     #   qs = 0.0 + 0.0im
     #   if θi < θh
     #       # Angle of the flat
     #       Δxf, Δyf = hills[2] - hills[1]
     #       θf = atan2(Δyf,Δxf)
     #       # Angle of the image path relative to the flat
     #       θ = θi - θf
     #       # Net image source to receiver height, as needed by QQ
     #       Δz = rr * sin(θ)
     #   
     #       rr = abs( rr * cos(θ) )
     #       qs = qq( rr, Δz, waveno, dbs[1] )
    #   end
    rr_Δz = calc_mirror( src_loc[2], hills, source=true )
    qs = isnothing(rr_Δz) ? 0.0 + 0.0im : qq( rr_Δz..., waveno, dbs[1]  )

 # Mirror Receiver
    #   # Distance from image receiver to top of hill
     #   Δx = rec_loc[2][1] - hills[3][1]
     #   Δy = hills[3][2] - rec_loc[2][2]
     #   rr = √( Δx^2 + Δy^2 )
     #   # Angle from image to top of hill
     #   θi = atan( Δy, Δx )
     #   # Angle of the hillside
     #   Δxh = hills[4][1] - hills[3][1]
     #   Δyh = hills[3][2] - hills[4][2]
     #   θh = atan( Δyh, Δxh )#
     #   
     #   qr = 0.0 + 0.0im
     #   if θi < θh
     #       # Angle of the flat
     #       Δxf = hills[5][1] - hills[4][1]
     #       Δyf = hills[4][2] - hills[5][2]
     #       θf = atan( Δyf, Δxf )
     #       # Angle of the image path relative to the flat
     #       θ = θi - θf
     #       # Net image source to receiver height, as needed by QQ
     #       Δz = rr * sin(θ)
     #       rr = abs( rr * cos(θ) )
     #       qr = ( rr, Δz, waveno, dbs[1] )
    #   end
    rr_Δz = calc_mirror( rec_loc[2], hills, source=false )
    qr = isnothing(rr_Δz) ? 0.0 + 0.0im : qq( rr_Δz..., waveno, dbs[1]  )

 # Wedge Angle
    Δx, Δz = hills[3] - hills[2]
    θ1 = atan( Δx, Δz )
    Δx = hills[4][1] - hills[3][1]
    Δz = hills[3][2] - hills[4][2]
    θ2 = atan( Δx, Δz )
    θ = θ1 + θ2

    pt = 0.0 + 0.0im
    f0dir = f1dir = f0refl = f1refl = 0
    for i in 1:4
        src_i = i == 1 || i == 3 ? 1 : 2
        rec_i = i <= 2 ? 1 : 2

        relev_refl_factor = i == 2 ? qs :
                                i == 3 ? qr :
                                    i == 4 ? qs * qr : 1.0 + 0.0im

        # Get length and angle of path from source to top of wedge
        Δx, Δz = hills[3] - src_loc[src_i]
        # Distance
        rh0 = √( Δx^2 + Δz^2 )
        # Angle, clockwise from straight down
        θh = atan( Δx, Δz )
        f0 = θh - θ1
        
        Δx = rec_loc[rec_i][1] - hills[3][1]
        Δz = hills[3][2] - rec_loc[rec_i][2]
        rh1 = √( Δx^2 + Δz^2 )
        # Angle, counterclockwise from straight down
        θ = atan( Δx, Δz )
        f1 = 2.0 * π - θ1 - θh

        if i == 1
            f0dir = f0
            f1dir = f1
        end
        if i == 4
            f0refl = f0
            f1refl = f1
        end

        tot_propag_path = rh0 + rh1

        if f0 > 0 && (f1 + θ) < (2.0 * π)
            h_over_wedgeleg = sin(f0) * tot_propag_path
            wedge_impedence1 = qq( tot_propag_path, h_over_wedgeleg, waveno, dbs[2] )
            h_over_wedgeleg = sin( 2 * π - f1 - θ ) * tot_propag_path
            wedge_impedence2 = qq( tot_propag_path, h_over_wedgeleg, waveno, dbs[2] )

            a = rh0 * rh1 / tot_propag_path
            any = 2.0 - θ / π
            # NON SONO CERTO I DUE MODI SIANO EQUIVALENTI
            #   pl = diffraction( tot_propag_path, a, f1-f0, -1.0, any, waveno ) +
            #        diffraction( tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1 +
            #        diffraction( tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2 +
            #        diffraction( tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2
            pl = sum( diffraction.(
                        tot_propag_path,
                        a,
                        [ f1-f0, f1+f0, f1+f0, f1-f0 ],
                        [  -1.0,  -1.0,   1.0,   1.0 ],
                        any,
                        waveno
                      ) .* [ 1.0, wedge_impedence1, wedge_impedence2, wedge_impedence1*wedge_impedence2 ] )
            pl *= relev_refl_factor
            pt += pl
        end
    end

 # Direct path source to receiver
    if π + f0dir - f1dir > 0
        Δx, Δz = rec_loc[1] - src_loc[1]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) / complex(rd, 0.0)
        pt += po
    end
 # Path mirrored source to receiver
    if π + f0refl - f1refl > 0
        Δx, Δz = rec_loc[1] - src_loc[2]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) * qq( rd, Δz, waveno, dbs[1] ) / complex(rd, 0.0)
        pt += po
    end   
 # Path source to mirrored receiver
    if π + f0refl - f1dir > 0
        Δx, Δz = rec_loc[2] - src_loc[1]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) * qq( rd, Δz, waveno, dbs[1] ) / complex(rd, 0.0)
        pt += po
    end

    Δx, Δz = rec_loc[1] - src_loc[1]
    rd = √( Δx^2 + Δz^2 )
    level = 4.34 * log( (rd * abs(pt))^2 )

    return level
end

function dal( first_second_dist::Real, second_third_dist::Real, src_h::Real, rec_h::Real, src_solpe_α::Real, flow_res1::Real, flow_res2::Real, freq::Real )::Real
    
    # Delany-Bazley for source and receiver leg
    dbs = delbaz.( freq, [flow_res1, flow_res2] )
    waveno = 2.0 * π * freq / 340.0

    calc_r( ra, rb, α ) = √( ra^2 + rb^2 - 2.0 * ra * rb * cos(α) )


    rh0 = √( src_h^2 + first_second_dist^2 )
    rh1 = √( rec_h^2 + second_third_dist^2 )
    r1 = rh0 + rh1
    f0 = atan(src_h, first_second_dist)
    f1 = 2.0 * π - src_solpe_α - atan(rec_h, second_third_dist)

    rd = calc_r( rh0, rh1, f0-f1 )
    direct_field = ℯ^complex(0.0, waveno*rd) / complex(rd, 0.0)

    h_over_wedgeleg = sin(f0) * r1
    wedge_impedence1 = qq( r1, h_over_wedgeleg, waveno, dbs[1] )
 #NON SO SE E' UN ERRORE
    # h_over_wedgeleg = sin(2.0 * π - f1 - src_solpe_α) * r1
    h_over_wedgeleg = sin(f1) * r1
    wedge_impedence2 = qq( r1, h_over_wedgeleg, waveno, dbs[2] )

    a = rh0 * rh1 / r1
    any = 2.0 - src_solpe_α / π
 # NON SONO CERTO I DUE MODI SIANO EQUIVALENTI
    #   pl = diffraction( tot_propag_path, a, f1-f0, -1.0, any, waveno ) +
    #        diffraction( tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1 +
    #        diffraction( tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2 +
    #        diffraction( tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2
    pl = sum( diffraction.(
                r1,
                a,
                [ f1-f0, f1+f0, f1+f0, f1-f0 ],
                [  -1.0,  -1.0,   1.0,   1.0 ],
                any,
                waveno
              ) .* [ 1.0, wedge_impedence1, wedge_impedence2, wedge_impedence1*wedge_impedence2 ] )

    if f1-f0 < π
        pl += direct_field
    end

    if f1+f0 < π
        rs = calc_r( rh0, rh1, f1+f0 )
        q = qq( rs, src_h+rh1*sin(π-f1), waveno, dbs[1] )
        ps = ℯ^complex(0.0, waveno*rs) / complex(rs, 0.0) * q
        pl += ps
    end

    θ = atan(rec_h, first_second_dist)

    if ( f1 - f0 + 2.0 * θ ) < π
        rr = clac_r( rh0, rh1, f1-f0+2.0*θ )
        q = qq(rr, rec_h+rh0*sin(f0 + src_solpe_α - π), waveno, dbs[2] )
        pr = ℯ^complex(0.0, waveno*rr) / complex(rr, 0.0) * q
        pl += pr
    end

    if ( f1 + f0 + 2.0 * θ ) < π
        rb = clac_r( rh0, rh1, f1+f0+2.0*θ )
        q1 = qq(rb, src_h+rh1*sin(2.0 * src_solpe_α - 3.0 * π + f1), waveno, dbs[1] )
        q2 = qq(rb, src_h+rh0*sin(-f0 + src_solpe_α - π), waveno, dbs[2] )
        pb = ℯ^complex(0.0, waveno*rb) / complex(rb, 0.0) * q1 * q2
        pl += pb
    end

    return 4.34 * log( (rd*abs(pl))^2 )
end

# FORSE VERSIONE MIGLIORATA
#=
function dal( first_second_dist::Real, second_third_dist::Real, src_h::Real, rec_h::Real, src_solpe_α::Real, flow_res1::Real, flow_res2::Real, freq::Real )::Real
    
    # Delany-Bazley for source and receiver leg
    dbs = delbaz.( freq, [flow_res1, flow_res2] )
    waveno = 2.0 * π * freq / 340.0

    rh0 = √( src_h^2 + first_second_dist^2 )
    rh1 = √( rec_h^2 + second_third_dist^2 )
    r1 = rh0 + rh1
    f0 = atan(src_h, first_second_dist)
    f1 = 2.0 * π - src_solpe_α - atan(rec_h, second_third_dist)
    θ = atan(rec_h, first_second_dist)

    h_over_wedgeleg = sin(f0) * r1
    wedge_impedence1 = qq( r1, h_over_wedgeleg, waveno, dbs[1] )
 #NON SO SE E' UN ERRORE
    # h_over_wedgeleg = sin(2.0 * π - f1 - src_solpe_α) * r1
    h_over_wedgeleg = sin(f1) * r1
    wedge_impedence2 = qq( r1, h_over_wedgeleg, waveno, dbs[2] )

    a = rh0 * rh1 / r1
    any = 2.0 - src_solpe_α / π
 # NON SONO CERTO I DUE MODI SIANO EQUIVALENTI
    #   pl = diffraction( tot_propag_path, a, f1-f0, -1.0, any, waveno ) +
    #        diffraction( tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1 +
    #        diffraction( tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2 +
    #        diffraction( tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2
    pl = sum( diffraction.(
                r1,
                a,
                [ f1-f0, f1+f0, f1+f0, f1-f0 ],
                [  -1.0,  -1.0,   1.0,   1.0 ],
                any,
                waveno
              ) .* [ 1.0, wedge_impedence1, wedge_impedence2, wedge_impedence1*wedge_impedence2 ] )

    αs = [ f1-f0, f1+f0, f1-f0+2θ, f1+f0+2θ ]
    βs = [ π-f1, f0+src_solpe_α-π, 2src_solpe_α-3π+f1, -f0+src_solpe_α-π ]
    dists = [ src_h+rh1, rec_h+rh0 ]
    for i in 1:4
        if αs[i] < π
            r = √( rh0^2 + rh1^2 - 2.0 * rh0 * rh1 * cos(αs[i]) )
            if i == 1
                q = 1
            elseif i == 4
                q = qq( r, dists[1]*sin(βs[i-1]), waveno*rα, db[1] ) * qq( r, dists[2]*sin(βs[i]), waveno*rα, db[2] )
            else
                q = qq( r, dists[2-i%2]*sin(βs[i-1]), waveno*rα, db[2] )
            end
            p = ℯ^complex(0.0, waveno*rα) / complex(rα, 0.0) * q
            pl += p
        end
    end

    return 4.34log( (rd*abs(pl))^2 )
end
=#



atten = onCut( distances_profiles[1], heights_profiles[1], zeros(length(heights_profiles)), dtm[r0,c0][1], heights_profiles[1][end], 1, [dB] )

atten = onCut( distances_profiles[5], heights_profiles[5], zeros(length(heights_profiles)), dtm[r0,c0][1], heights_profiles[1][end], 1, [dB] )

atten = onCut( distances_profiles[600], heights_profiles[600], zeros(length(heights_profiles)), dtm[r0,c0][1], heights_profiles[1][end], 1, [dB] )

atten = onCut( distances_profiles[201], heights_profiles[201], zeros(length(heights_profiles)), dtm[r0,c0][1], heights_profiles[1][end], 1, [dB] )



function onCut( distances::AbstractVector, heights::AbstractVector, impdcs::AbstractVector, src_h::Real, rec_h::Real, nfreq::Int64, freqs::AbstractVector )

    ihard = 0
    isoft = 0 
    flow_ress = [0.0, 0.0]
    atten = 0.0
    attenuations = []
    profile = []
    flowpr = []
    ignd = []


 # ====================================================== Section to Process Profile =====================================================================

    for i in 1:length(distances)
        push!( profile, (distances[i], heights[i]) )
        push!( flowpr, (distances[i], impdcs[i]) )
        if flowpr[i][2] <= 1000.0
            push!( ignd, 0 )
            isoft += 1
            flow_ress[1] += flowpr[i][2]
        else
            push!( ignd, 1 )
            ihard += 1
            flow_ress[2] += flowpr[i][2]
        end
    end
    flow_ress[1] /= max(isoft, 1)
    flow_ress[2] /= max(ihard, 1)
    
    if flow_ress[1] == 0.0
        flow_ress[1] = 200.0
    end
    if flow_ress[2] == 0.0
        flow_ress[2] = 10.0^6
    end

 # =======================================================================================================================================================

    # Points of interest:
    #     min near src ||   max  || min near rec
    #           x   z  || x   z  || x   z
    points = [ (0, 0.0), (0, 0.0), (0, 0.0) ]
    # Max point (second point)
    points[2] = minmax( profile, src_h, rec_h )[2]
    # Min point near source (first point)
    points[1] = minmax( profile[ 1:points[2][1] ], src_h, 0.0 )[1]
    # Min point near receiver (third point)
    ntemp = length(profile) - points[2][1]
    #   xs = [ point[1] for point in profile[ points[2][1]+1:ntemp ] ] 
    #   res = minmax( xs, 0.0, rec_h )[1]
    res = minmax( profile[ points[2][1]+1:ntemp ], 0.0, rec_h )[1]
    points[3] = (points[2][1]-1, 0.0) + res

    hillxz = [ profile[1] ]
    klocs = [1]
    nmm = 1
    for i in 1:3
        push!( hillxz, profile[ points[i][1] ] )
 # NEL CODICE METTONO CHE SE "points[i][1] == klocs[nmm]" VANNO AL CONTINUE
  # SCRITTO ALLA FINE DEL CICLO, NON SO SE VOGLIA DIRE CHE IN QUEL CASO SALTA
  # TUTTO IL CICLO
        if points[i][1] != klocs[nmm]
            nmm += 1
            if length(klocs) >= nmm
                klocs[nmm]= points[i][1]
            else
                push!( klocs, points[i][1] )
            end
        end
    end

    push!( hillxz, profile[end] )
    if klocs[nmm] < length(profile)
        nmm += 1
        if length(klocs) >= nmm
            klocs[nmm] = length(profile)
        else
            push!( klocs, length(profile) )
        end
    end

    if points[1][1] == points[2][1]
        hillxz[2] = hillxz[1]
    end

    if points[3][1] == points[2][1]
        hillxz[4] = hillxz[5]
    end

    if hillxz[1][1] == hillxz[2][1]
        dx, dy = hillxz[3] - hillxz[1]
        hillxz[2] = hillxz[1] + 0.1(dx, dy)
    end

    if hillxz[4][1] == hillxz[5][1]
        dx, dy = hillxz[5] - hillxz[3]
        hillxz[4] = hillxz[3] + 0.1(dx, dy)
    end

 # ================================================= Profile and Supporting Stuff Established ============================================================
 
  # ----------------------------------------------- Level Model ----------------------------------------------------- 
    
    if nmm <= 2
        println("nmm <= 2: $(nmm <= 2)" )
        ax = profile[1][1]
        ay = 0.0
        az = profile[1][2]
        ox = profile[end][1]
        oy = 0.0
        oz = profile[end][2]

        d = √( (ax-ox)^2 + (ay-oy)^2 )

        for j in 1:nfreq
            duml = 0.0
            if ihard > 0
    # NON MI E' CHIARO IL SENSO DI "duml"
                duml, atten = egal( d/2, d/2, src_h, rec_h, flow_ress[2], flow_ress[2], 0.0, 0.0, 0.0, freq[j], duml )
                atten = attenh
            end
            if isoft > 0
    # NON MI E' CHIARO IL SENSO DI "duml"
                duml, attens = egal( d/2, d/2, src_h, rec_h, flow_ress[1], flow_ress[1], 0.0, 0.0, 0.0, freq[j], duml )
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

  # ------------------------------------------------- Hill Model ---------------------------------------------------- 
    
    zz2 = 0
    zcrit = 0
    if nmm == 3
        println("nmm == 3: $(nmm == 3)" )
        # Total distance dist and distance to second point dsl
        dist = profile[ klocs[3] ][1] - profile[ klocs[1] ][1]
        dsl = profile[ klocs[2] ][1] - profile[ klocs[1] ][1]
        # Altitudes at the three points
        zz1 = profile[ klocs[1] ][2]
        zz2 = profile[ klocs[2] ][2]
        zz3 = profile[ klocs[3] ][2]
        zcrit = zz1 + (zz3-zz1) * dsl / dist
    end
    
    if nmm > 3 || ( nmm == 3 && zz2 >= zcrit )
        println("nmm > 3 || ( nmm == 3 && zz2 >= zcrit ): $(nmm > 3 || ( nmm == 3 && zz2 >= zcrit ))" )
        # Set up the source and receiver locations, normal to the corresponding plateaus
        cosθ = 1.0
        sinθ = 0.0
        Δx, Δz = hillxz[2] - hillxz[1]
        if Δx > 0
            θ = atan(Δz, Δx)
            cosθ = cos(θ)
            sinθ = sin(θ)
        end

        # Source location is hs above the start of the terrain cut
        srcloc = [ hillxz[1] + (0, src_h) ]
        # Reflect the original source image
        push!( srcloc, srcloc[1] + 2*src_h*cosθ*(sinθ, -cosθ) )

        cosθ = 1.0
        sinθ = 0.0
        Δx, Δz = hillxz[5] - hillxz[4]
        if Δx > 0
            θ = atan(Δz, Δx)
            cosθ = cos(θ)
            sinθ = sin(θ)
        end

        # Right over the end of the cut receiver
        recloc = [ hillxz[5] + (0, rec_h) ]
        # reflect the original receiver image
        push!( recloc, recloc[1] + 2*rec_h*cosθ*(sinθ, -cosθ) )

        for j in 1:nfreq
            if ihard > 0
                attenh = bakkernn( srcloc, recloc, hillxz, flow_ress[2], flow_ress[2], flow_ress[2], freqs[j] )
                atten = attenh
                println("Atten H: $attenh")
            end
            if isoft > 0
                attens = bakkernn( srcloc, recloc, hillxz, flow_ress[1], flow_ress[1], flow_ress[1], freqs[j] )
                atten = attens
                println("Atten S: $attens")
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

  # ------------------------------------------------ Valley Model --------------------------------------------------- 

    if nmm == 3 && zz2 < zcrit
        println("nmm == 3 && zz2 < zcrit: $(nmm == 3 && zz2 < zcrit)" )
        diff1 = profile[ klocs[2] ] - profile[ klocs[1] ]
        diff2 = profile[ klocs[3] ] - profile[ klocs[2] ]

        α1 = atan( diff1[2], diff1[1] )
        α2 = atan( diff2[2], diff2[1] )
        α = α2 - α1 + π
        # d0 = √( diff1[1]^2 + diff1[2]^2 )
        # d1 = √( diff2[1]^2 + diff2[2]^2 )
        d0, d1 = @. √sum( [diff1, diff2]^2 )

        for j in 1:nfreq
            attenh = 0.0
            attens = 0.0
            if ihard > 0
                attenh = dal( d0, d1, src_h, rec_h, α, flow_ress[2], flow_ress[2], freqs[j] )
                atten = attenh
                println("Atten H: $attenh")
            end
            if isoft > 0
                attenh = dal( d0, d1, src_h, rec_h, α, flow_ress[1], flow_ress[1], freqs[j] )
                atten = attens
                println("Atten S: $attens")
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

    println("Nmm: $nmm")
    return attenuations
end







# Taken from "https://www.geeksforgeeks.org/dda-line-generation-algorithm-computer-graphics/"
function DDA( map, x0::Number, y0::Number, xn::Number, yn::Number )
    Δx = xn - x0
    Δy = yn - y0
    steps = max( abs(Δx), abs(Δy) )
    x_inc = Δx / steps
    y_inc = Δy / steps
    x = x0
    y = y0
    heigths_profile = []
    coords_profile = []
    for i in 1:steps
        push!( heigths_profile, map[toInt(x), toInt(y)][1] )
        push!( coords_profile, Tuple(GeoArrays.coords( map, [toInt(x), toInt(y)])) )
        x += x_inc
        y += y_inc
    end
    return heigths_profile, coords_profile
end



#   https://tildesites.bowdoin.edu/~ltoma/teaching/cs350/spring06/Lecture-Handouts/gis-viewshedsKreveld.pdf
#=
    function ground_loss( x0::Real, y0::Real, dB::Real, heights_map, impedences_map )
    end
    # Punto 3870, 4420 (centro):
    #   x = 723204.0
    #   y = 5065493.0
    #   ground_loss( x, y, 110, *( @__DIR__, "\\Mappe\\DTM_32.tiff" ) )
=#

x0 = 710705.0
y0 = 5065493.0
dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
dtm = GeoArrays.read(dtm_file)
# Coordinate dei punti
x1, y1 = GeoArrays.coords( dtm, [1,1] )
xn, yn = GeoArrays.coords( dtm, size(dtm)[1:2] )
# Dimensioni in metri di una cella
Δx, Δy = (xn-x1, y1-yn) / size(dtm)[1:2]

dB = 110 - 32
r0, c0 = GeoArrays.indices( dtm, [x0, y0] )
max_radius = ceil(10^(dB/20))
cell_num = toInt( max_radius / Δx, ceil )

heights_profiles = [
    dtm[ r0, c0-cell_num:c0 ], # x axis 1
    dtm[ r0, c0:c0+cell_num ], # x axis 2
    dtm[ r0-cell_num:r0, c0 ], # y axis 1
    dtm[ r0:r0+cell_num, c0 ], # y axis 2
    # I VALORI OTTENUTI IN QUESTO MODO NON SONO IN ORDINE
    # MANCANO I VALORI DI x0y0
    [ dtm[r0+i, c0-i][1] for i in (-cell_num):0 ], # first diagonal 1
    [ dtm[r0+i, c0-i][1] for i in 0:cell_num ], # first diagonal 2
    [ dtm[r0+i, c0+i][1] for i in (-cell_num):0 ], # second diagonal 1
    [ dtm[r0+i, c0+i][1] for i in 0:cell_num ]  # second diagonal 2
    # vcat( [ dtm[r0+i, c0-i] for i in 1:cell_num ],  [ dtm[r0-i, c0+i] for i in 1:cell_num ] ), # first diagonal
    # vcat( [ dtm[r0-i, c0-i] for i in 1:cell_num ],  [ dtm[r0+i, c0+i] for i in 1:cell_num ] ) # second diagonal
]
coords_profiles = [
    [ GeoArrays.coords(dtm, [r0, col]) for col in c0-cell_num:c0 ],
    [ GeoArrays.coords(dtm, [r0, col]) for col in c0:c0+cell_num ],
    [ GeoArrays.coords(dtm, [row, c0]) for row in r0-cell_num:r0 ],
    [ GeoArrays.coords(dtm, [row, c0]) for row in r0:r0+cell_num ],
    [ GeoArrays.coords(dtm, [r0+i, c0-i]) for i in (-cell_num):0 ],
    [ GeoArrays.coords(dtm, [r0+i, c0-i]) for i in 0:cell_num ],
    [ GeoArrays.coords(dtm, [r0+i, c0+i]) for i in (-cell_num):0 ],
    [ GeoArrays.coords(dtm, [r0+i, c0+i]) for i in 0:cell_num ],
]

row_begin = r0 - cell_num
row_end = r0 + cell_num
col_begin = c0 - cell_num
col_end = c0 + cell_num


#= CREAZIONE PROFILI SU DUE QUADRANTI
    # Upper side of the square, plus the halves of left and right side above x axis
    top_idxs = [ (row_begin, col) for col in  col_begin:col_end ]
    left_idxs = [ (row, col_begin) for row in row_begin+1:r0-1 ]
    right_idxs = [ (row, col_end) for row in  row_begin+1:r0-1 ]

    # Vector with the indexis of the borders of the area of effect above the x axis
    endpoints = vcat( right_idxs, top_idxs, left_idxs )


    # Rotazioni:
     #  xa₂ = x0 + (xa - x0)cos(θ) - (ya - y0)sin(θ)
     #  ya₂ = y0 + (xa - x0)sin(θ) - (ya - y0)cos(θ)

    # TUTTI GLI ARRAY TRANNE QUELLI MESSI MANUALMENTE HANNO UN VALORE IN MENO
    # For each point compute the rasterization of the line that connects it to its symmetrical point using the center (r0 c0) as reference
    for point in endpoints
        # If the point is on the y axis or on the diagonals skip ( the points cannot be on the x axis by construction of the side indexis vectors )
        if point[2] == c0 || abs(point[1] - r0) == abs(point[2] - c0)
            continue
        end
        # Compute profile from center to point
        heights, coords = DDA( dtm, r0, c0, point[1], point[2] )
        push!( heights_profiles, heights )
        push!( coords_profiles, coords )
        # Compute symmetrical profile
        rm, cm = 2(r0, c0) - point
        heights, coords = DDA( dtm, r0, c0, rm, cm )
        push!( heights_profiles, heights )
        push!( coords_profiles, coords )
    end
=#


top_idxs = [ (row_begin, col) for col in  c0+1:col_end ]
right_idxs = [ (row, col_end) for row in  row_begin+1:r0-1 ]
# Vector containing the indexis of the points on first quadrant of the border of the area of interest
endpoints = vcat( top_idxs, right_idxs )
# For every aforementioned points compute the profile ranging from them and their rotations to the center point
for point in endpoints
    if point[2] == c0 || abs(point[1] - r0) == abs(point[2] - c0)
        continue
    end
    # Compute profile from center to point
    heights, coords = DDA( dtm, r0, c0, point[1], point[2] )
    push!( heights_profiles, heights )
    push!( coords_profiles, coords )
    # Compute the profiles for the three points obtained by rotating the starting point by 90, 180 and 270 degrees respectively
    for α in [90, 180, 270]
        rm, cm = rotate_point( point[1], point[2], r0, c0, α )
        heights, coords = DDA( dtm, r0, c0, rm, cm )
        push!( heights_profiles, heights )
        push!( coords_profiles, coords )
    end
end

distances_profiles = map.( point -> √((point[1] - x0)^2 + (point[2] - y0)^2), coords_profiles )

atten = onCut( distances_profiles[1], heights_profiles[1], zeros(length(heights_profiles)), dtm[r0,c0][1], heights_profiles[1][end], 1, [dB] )











attenuations = []
for i in 1:length(heights_profiles)
    atten = onCut( distances_profiles[i], heights_profiles[i], zeros(length(heights_profiles)), dtm[r0,c0][1], heights_profiles[1][end], 1, [dB] )
    push!( attenuations, atten )
end









#= COSE CON GEODATA
    using GeoData

    #   dtm_file = *( @__DIR__, "\\Mappe\\DTM_wgs84.tiff" )
    dtm_file = *( @__DIR__, "\\Mappe\\DTM_32.tiff" )
    dtm = GeoArray(dtm_file)
    plot(dtm)
    #   # Altezza di x, y:
     #   dtm[x, y]
     #   # Coordinate:
     #    # X (minima e massima)
     #   dtm.dims[1][1]
     #   dtm.dims[1][end]
     #   length(dtm.dims[1])
     #    # Y (massima e minima)
     #   dtm.dims[2][1]
     #   dtm.dims[2][end]
    #   length(dtm.dims[2])

    Δs = ( dtm.dims[1][end], dtm.dims[2][1] )
    Δs -= ( dtm.dims[1][1], dtm.dims[2][end] )
    Δs /= Float64.(size(dtm)[1:2])
    dB = 110 - 32
    max_radius = ceil(10^(dB/20))
=#

#=  CREAZIONE GRIGLIA
    function create_grid( xs::Real, ys::Real, dB::Real )
        # maximum ideal radius for the diffusion of a sound with given intensity  
        max_radius = ceil(10^(dB/20))
        # Coordinates of the first point (upper left) of the grid containing the radius 
        x0 = xs - max_radius 
        y0 = ys + max_radius

        grid = CartesianGrid( (x0, ys-max_radius), (xs+max_radius, y0), (200.0,200.0) )  

    #=
        # LE COORDINATE DELL'AREA NON SONO VALIDE PERCHE' NON SONO ADEGUATAMENTE SCALATE
        reg = RectRegion( "NEW", "VNT", "Area Of Interest", [ x0, ys-max_radius, xs+max_radius, y0 ] )

        return grid, reg
    =#
        return grid
    end

    g = create_grid( 22., 22., 56. )
    plot(g)
    centers = centroid.(g)
    plot!(centers)

    size = convert( Int64, √length(centers) )
    start = size * Int64(floor(size/2)) + 1
    x_axis = centers[ start:start+size-1 ]
    plot!(x_axis, c=:blue )

    start = convert( Int64, ceil(size/2) )
    y_axis = centers[[ start + size*i for i in 0:size-1 ]]
    plot!(y_axis, c=:red )

    diagonal1 = centers[[ size + (size-1)i for i in 0:size-1 ]]
    plot!(diagonal1, c=:yellow)

    diagonal2 = centers[[ 1 + (size+1)i for i in 0:size-1 ]]
    plot!(diagonal2, c=:purple)

    is = collect(x0:size)
    x0 = y0 = Int64(round(size/2))
    r = Int64(floor(size/2))
    radius_y(x) = Int64( round( √( r^2 - (x-x0)^2 ) + y0 ) )
    idxs = ( radius_y.(is) .- 1 )
    half_circle = centers[idxs]
    plot!( half_circle, c=:orange )
=#

#=  PLOTTING DI DTM E ALTRO
    using ArchGDAL
    using GeoStats
    using GeoRegions
    using Plots
    using Shapefile

    #   geo = GeoRegion("GLB")  # Rappresenta il mondo
    #   bound_lon, bound_lat = coordGeoRegion(geo)  # Ottieni i vettori di latitudine e longitudine se `geo` e poligonale ritorna anche i vettori di lat e lon per la forma
    #   ref = ( max(bound_lon...) + min(bound_lon...) ) / 2
    #   ginfo = RegionGrid(geo, bound_lon, bound_lat )
    #   plot!( bound_lon.-ref, bound_lat )

    # Shapefile
    coast_file = *( @__DIR__, "\\Mappe\\ne_10m_coastline\\ne_10m_coastline.shp" )
    shp_tb = Shapefile.Table(coast_file)
    shp = Shapefile.shapes(shp_tb)
    plot!(shp)

    # Shapefile Veneto
    veneto_shp = *( @__DIR__, "\\Mappe\\c0104011_Comuni\\c0104011_Comuni.shp")
    shp_tb2 = Shapefile.Table(veneto_shp)
    shp2 = Shapefile.shapes(shp_tb2)
    plot!(shp2)

    # Poligono veneto
    veneto = try
                GeoRegion("VNT")
             catch
                # [ 9.5 47.0, 14.0 47.0, 14.0 44.0, 9.5 44.0, 9.5 47.0 ]
                RectRegion( "VNT", "GLB", "Veneto", [ 47.0, 44.0, 14.0, 9.5 ] )
             end
    blon, blat = coordGeoRegion(veneto)
    plot!( blon, blat )


    # LA GRIGLIA E' FUORI SCALA RISPETTO A TUTTO IL RESTO
    geo_grid, aoi = create_Grid( 22., 22., 56. )
    plot!( geo_grid )



    plot!( [ 22.0, 22.0 ] )




    # Tiff file 
    dtm_file = *( @__DIR__, "\\Mappe\\DTM_32.tif" )
    map_file = *( @__DIR__, "\\Mappe\\ne_10m_coastline" )
    map_v_file = *( @__DIR__, "\\Mappe\\c0104011_Comuni" )

    dtm = ArchGDAL.read(dtm_file)
    map = ArchGDAL.read(map_file)
    map_v = ArchGDAL.read(map_v_file)

    band = ArchGDAL.getband( dtm, 1 )
    ArchGDAL.imread(band)
    ArchGDAL.imread(map)
    ArchGDAL.imread(map_v)
    plot(map)
=#

end # module