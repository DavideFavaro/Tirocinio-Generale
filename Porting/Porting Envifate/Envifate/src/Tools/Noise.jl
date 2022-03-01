module Noise
"""
"""



using ArchGDAL
using GeoArrays
using Plots
using Shapefile


include("..\\Library\\Functions.jl")


const agd = ArchGDAL
const ga = GeoArrays



repeat!(A::AbstractVector, count::Integer ) = append!( A, repeat(A, count-1) )
# Functions to find the coordinates of the point resulting from the rotation of "(xp, yp)" by a angle "θ" around "(xc, yc)"
rotate_x( xp, yp, xc, yc, θ ) = round( Int64, (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) + xc )
rotate_y( xp, yp, xc, yc, θ ) = round( Int64, (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) + yc )
rotate_point( xp, yp, xc, yc, θ ) = θ == 0 ? (xp, yp) : ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )
distance( p0::Tuple{Float64, Float64}, pn::Tuple{Float64, Float64} ) =  √( sum( v -> v^2, pn .- p0 ) )



"""
    transmission_loss( r::Real )

Compute the transmission loss of a noise over `r` distance
"""
function transmission_loss( radius::Real )
    return 20log10(radius)
end

function atmospheric_absorpion_loss( radius::Real, height_m::Real, relative_humidity::Real, temperature_k::Real, frequency::Real )
    # Calculate atmospheric absorption coefficient using ANSI S1.26-1995 standard
    # Convert elevation to atmospheric pressure
    atmospheric_pressure = 101.325( 1 - ( 2.25577 * 10^(-5) * height_m ) )^5.25588
    # Calculate derived values for subsequent equations
    p_atm_pressure = atmospheric_pressure / 101.325
    t_tr = temperature_k / 293.15
    # Convert relative humidity to molar concentration of water vapor
    C = ( -6.8346( 273.16 / temperature_k )^1.261 ) + 4.6151
    p_saturation_pressure = 10^C
    humidity = relative_humidity * p_saturation_pressure * p_atm_pressure^(-1)
    # Calculate relaxation frequency of O (equation 3)
    #   frO₂ = ( p_atm_pressure * ( (24 + 4.04e04) * humidity ) * (0.02 + humidity) ) / (0.391 + humidity)
    frO₂ = p_atm_pressure * ( 24 + ( 4.04e04humidity * ( (0.02 + humidity) / (0.391 + humidity) ) ) )
    # Calculate relaxation frequency of N (equation 4)
    frN₂ = p_atm_pressure * √t_tr * ( 9 + ( 280humidity * ℯ^( -4.170(t_tr^(-1/3) - 1) ) ) )
    # Calculate alpha (equation 5)
    term1 = 1.84 * 10^(-11) * p_atm_pressure^(-1) * √t_tr
    #   term2 = t_tr^(-2.5) * ( 0.01275 * ℯ^(-2239.1 / temp_k) * ( frO₂ / (frO₂^2 + freq^2) ) )
    term2 = t_tr^(-2.5)( 0.01275ℯ^(-2239.1 / temp_k) / ( frO₂ + (freq^2 / frO₂) ) )
    #   term3 = 0.1068 * ℯ^(-3352 / temp_k) * ( frN₂ / (frN₂^2 + freq^2) )
    term3 = t_tr^(-2.5) * ( 0.1068ℯ^(-3352 / temp_k) / ( frN₂ + (freq^2 / frN₂) ) ) 
    #   α = 8.686 * (frequency^2) * ( term1 + term2 + term3 )
    α = frequency^2 * (term1 + term2 + term3)

    return α * radius / 100
end

function minmax( profile::Vector, rel_h_src::Real=0.0, rel_h_rec::Real=0.0 )::Vector{Tuple{Int64, Real}}
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
        m = round( transition_height / ( la / 6.0 ) )
        fixed_speed = 2.0 * π * freq / (
                                 340.0 + (
                                             ( m-1 * la/6.0 + la/10.0 ) / 2.0 
                                         ) * e_wind_vel/10.0
                             )
        e_turbulence_scale = 10e3 / freq

        eq(x, rrn, rrm ) = ℯ^complex( 0.0, -x*(rrn + rrm) )
        eq_cr( rrn, rrm ) = rrm * √( rrm * rrn * (rrn + rrm) )
        eq_jarr( kv, k, rrn, rrm; q1=nothing, q2=nothing ) =  ( eq(kv, rrn, rrm) - eq(k, rrn, rrm) ) * ( isnothing(q1) ? 1 : q1 ) * ( isnothing(q2) ? 1 : q2 ) / complex( eq_cr(rrn, rrm), 0.0 ) 

        for j in 1:m
            ha = (j-1) * (la/6.0) + (la/10.0)
            v = [ e_wind_vel + turbulence / 10.0 * cos( ha * 2π/e_turbulence_scale ) ]
            push!( 
                v,
                2e_wind_vel - v[1],
                e_wind_vel
            )

            for i in 1:3
                kv[i] = 2.0 * π * freq / ( 340.0 + (ha/2.0) * v[i]/10.0 )

                rr = [
                    √( d1^2 + (src_h - ha)^2 ),
                    √( d1^2 + (src_h + ha)^2 ),
                    √( d2^2 + (ha - rec_h)^2 ),
                    √( d2^2 + (ha + rec_h)^2 ) 
                ]

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

    end
    q2 = qq2( d1+d2, src_h, rec_h, freq, rec_rx )
    c = ℯ^complex( 0.0, -fixed_speed*rr ) / complex( rr, 0.0 ) * q2
    ch = ℯ^complex( 0.0, -fixed_speed*rd ) / complex( rd, 0.0 )
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

#= VERSIONE CON IL CALCOLO DELLE ATTENUAZIONI PARZIALI (FORSE)
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

    results = []
    Δl = sum = denom = 0
    for i in 1:length(dists)
        w1 = w2 = 0
        if dists[i] > srcInf && dists[i] < recInf
            w1 = dists[i] < drefl ? rec_h : src_h
            w2, Δl = ground_type[i] == Hard ? (0.5, hard_atten) : (1.0, soft_atten)
        end
        sum += w1 * w2 * Δl
        denom += w1 * w2
        push!( results, sum/denom )
    end

    return results
end
=#
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

function diffraction!( aalast::Ref{Float32}, r1::Real, a::Real, al2::Real, pm::Real, any::Real, k::Real )::Complex
    df = -ℯ^complex(0.0, k * r1 + π / 4.0) / complex(r1, 0.0)
    tangent = tan( (π + pm * al2) / (2.0 * any) )
    aa = tangent != 0 ? ( 1.0 / tangent / (2.0 * any) / √(2.0 * π * k * a) ) : aalast[]
    aalast[] = aa
    df *= complex(aa, 0.0)
    n = pm < -0.9 ?
            al2 > π + any * π ? 1 :
                al2 > π - any * π ? 0 : -1 :
            al2 > any * π - π ? 1 : 0
    xv = 2.0 * k * a * cos( (2.0 * n * any * π - al2) / 2.0 )^2
    return df * ℯ^complex(0.0, -xv) * fres(√xv) * complex(0.0, -2.0√xv)
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

function bakkernn( hills::AbstractVector, src_loc::AbstractVector, rec_loc::AbstractVector, src_flow_res::Real, ber_flow_res::Real, rec_flow_res::Real, freq::Real )::Real
    # Delany-Bazley under source, berm and reciver
    dbs = delbaz.( freq, [ src_flow_res, ber_flow_res, rec_flow_res ] )
    waveno = 2.0 * π * freq / 340.0

 # Mirror Source
    # Distance from image source to top of hill
    Δx, Δy = hills[3] .- src_loc[2]
    rr = √( Δx^2 + Δy^2 )
    # Angle from image to top of hill
    θi = atan( Δy, Δx )
    # Angle of the hillside
    Δxh, Δyh = hills[3] .- hills[2]
    θh = atan( Δyh, Δxh )
    
    qs = 0.0 + 0.0im
    if θi < θh
        # Angle of the flat
        Δxf, Δyf = hills[2] .- hills[1]
        θf = atan( Δyf, Δxf )
        # Angle of the image path relative to the flat
        θ = θi - θf
        # Net image source to receiver height, as needed by QQ
        Δz = rr * sin(θ)
    
        rr = abs( rr * cos(θ) )
        qs = qq( rr, Δz, waveno, dbs[1] )
    end
    #   rr_Δz = calc_mirror( src_loc[2], hills, source=true )
    #   qs = isnothing(rr_Δz) ? 0.0 + 0.0im : qq( rr_Δz..., waveno, dbs[1]  )

 # Mirror Receiver
    # Distance from image receiver to top of hill
    Δx = rec_loc[2][1] - hills[3][1]
    Δy = hills[3][2] - rec_loc[2][2]
    rr = √( Δx^2 + Δy^2 )
    # Angle from image to top of hill
    θi = atan( Δy, Δx )
    # Angle of the hillside
    Δxh = hills[4][1] - hills[3][1]
    Δyh = hills[3][2] - hills[4][2]
    θh = atan( Δyh, Δxh )
    
    qr = 0.0 + 0.0im
    if θi < θh
        # Angle of the flat
        Δxf = hills[5][1] - hills[4][1]
        Δyf = hills[4][2] - hills[5][2]
        θf = atan( Δyf, Δxf )
        # Angle of the image path relative to the flat
        θ = θi - θf
        # Net image source to receiver height, as needed by QQ
        Δz = rr * sin(θ)
        rr = abs( rr * cos(θ) )
        qr = qq( rr, Δz, waveno, dbs[1] )
    end
    #   rr_Δz = calc_mirror( rec_loc[2], hills, source=false )
    #   qr = isnothing(rr_Δz) ? 0.0 + 0.0im : qq( rr_Δz..., waveno, dbs[1]  )

 # Wedge Angle
    Δx, Δz = hills[3] .- hills[2]
    θ1 = atan( Δx, Δz )
    Δx = hills[4][1] - hills[3][1]
    Δz = hills[3][2] - hills[4][2]
    θ2 = atan( Δx, Δz )
    θ = θ1 + θ2

    pt = 0.0 + 0.0im
    f0dir = 0
    f1dir = 0
    f0refl = 0
    f1refl = 0
    for i in 1:4
        src_i = ( i == 1 || i == 3 ) ? 1 : 2
        rec_i = i <= 2 ? 1 : 2

        relev_refl_factor = i == 2 ? qs :
                                i == 3 ? qr :
                                    i == 4 ? qs * qr : 1.0 + 0.0im

        # Get length and angle of path from source to top of wedge
        Δx, Δz = hills[3] .- src_loc[src_i]
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
        f1 = 2π - θ1 - θh

        if i == 1
            f0dir = f0
            f1dir = f1
        end
        if i == 4
            f0refl = f0
            f1refl = f1
        end

        tot_propag_path = rh0 + rh1

        if f0 > 0 && ( (f1 + θ) < 2.0π )
            h_over_wedgeleg = sin(f0) * tot_propag_path
            wedge_impedence1 = qq( tot_propag_path, h_over_wedgeleg, waveno, dbs[2] )
            h_over_wedgeleg = sin( 2 * π - f1 - θ ) * tot_propag_path
            wedge_impedence2 = qq( tot_propag_path, h_over_wedgeleg, waveno, dbs[2] )

            a = rh0 * rh1 / tot_propag_path
            any = 2.0 - θ / π
            aalast = Ref{Float32}(0.0f0)
            dif1 = diffraction!( aalast, tot_propag_path, a, f1-f0, -1.0, any, waveno )
            dif2 = diffraction!( aalast, tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1
            dif3 = diffraction!( aalast, tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2
            dif4 = diffraction!( aalast, tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2
            pl = dif1 + dif2 + dif3 + dif4

            # NON SONO CERTO I DUE MODI SIANO EQUIVALENTI
            #   pl = diffraction( tot_propag_path, a, f1-f0, -1.0, any, waveno ) +
            #        diffraction( tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1 +
            #        diffraction( tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2 +
            #        diffraction( tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2

            #   pl = sum( diffraction.(
            #               tot_propag_path,
            #               a,
            #               [ f1-f0, f1+f0, f1+f0, f1-f0 ],
            #               [  -1.0,  -1.0,   1.0,   1.0 ],
            #               any,
            #               waveno
            #             ) .* [ 1.0, wedge_impedence1, wedge_impedence2, wedge_impedence1*wedge_impedence2 ] )
            pl *= relev_refl_factor
            pt += pl
        end
    end

 # Direct path source to receiver
    if π + f0dir - f1dir > 0
        Δx, Δz = rec_loc[1] .- src_loc[1]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) / complex(rd, 0.0)
        pt += po
    end
 # Path mirrored source to receiver
    if π + f0refl - f1dir > 0
        Δx, Δz = rec_loc[1] .- src_loc[2]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) * qq( rd, Δz, waveno, dbs[1] ) / complex(rd, 0.0)
        pt += po
    end   
 # Path source to mirrored receiver
    if π + f0dir - f1refl > 0
        Δx, Δz = rec_loc[2] .- src_loc[1]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) * qq( rd, Δz, waveno, dbs[1] ) / complex(rd, 0.0)
        pt += po
    end

    Δx, Δz = rec_loc[1] .- src_loc[1]
    rd = √( Δx^2 + Δz^2 )
    level = 4.34log( (rd * abs(pt))^2 )

    return level
end

function dal( first_second_dist::Real, second_third_dist::Real, src_h::Real, rec_h::Real, src_slope_α::Real, flow_res1::Real, flow_res2::Real, freq::Real )::Real
    # Delany-Bazley for source and receiver leg
    dbs = delbaz.( freq, [flow_res1, flow_res2] )
    waveno = 2.0 * π * freq / 340.0

    calc_r( ra, rb, α ) = √( ra^2 + rb^2 - 2.0 * ra * rb * cos(α) )

    rh0 = √( src_h^2 + first_second_dist^2 )
    rh1 = √( rec_h^2 + second_third_dist^2 )
    r1 = rh0 + rh1
    f0 = atan(src_h, first_second_dist)
    f1 = 2.0 * π - src_slope_α - atan(rec_h, second_third_dist)

    rd = calc_r( rh0, rh1, f0-f1 )
    direct_field = ℯ^complex(0.0, waveno*rd) / complex(rd, 0.0)

    h_over_wedgeleg = sin(f0) * r1
    wedge_impedence1 = qq( r1, h_over_wedgeleg, waveno, dbs[1] )
    h_over_wedgeleg = sin(2.0 * π - f1 - src_slope_α) * r1
    wedge_impedence2 = qq( r1, h_over_wedgeleg, waveno, dbs[2] )

    a = rh0 * rh1 / r1
    any = 2.0 - src_slope_α / π
    aalast = Ref{Float32}(0.0)
 # NON SONO CERTO I DUE MODI SIANO EQUIVALENTI
    #   pl = diffraction( tot_propag_path, a, f1-f0, -1.0, any, waveno ) +
    #        diffraction( tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1 +
    #        diffraction( tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2 +
    #        diffraction( tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2
    pl = sum( diffraction!.(
                Ref(aalast),
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
        rr = calc_r( rh0, rh1, f1-f0+2.0*θ )
        q = qq(rr, rec_h+rh0*sin(f0 + src_slope_α - π), waveno, dbs[2] )
        pr = ℯ^complex(0.0, waveno*rr) / complex(rr, 0.0) * q
        pl += pr
    end

    if ( f1 + f0 + 2.0 * θ ) < π
        rb = calc_r( rh0, rh1, f1+f0+2.0*θ )
        q1 = qq(rb, src_h+rh1*sin(2.0 * src_slope_α - 3.0 * π + f1), waveno, dbs[1] )
        q2 = qq(rb, src_h+rh0*sin(-f0 + src_slope_α - π), waveno, dbs[2] )
        pb = ℯ^complex(0.0, waveno*rb) / complex(rb, 0.0) * q1 * q2
        pl += pb
    end

    return 4.34log( (rd*abs(pl))^2 )
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

function onCut( distances::Union{Vector{T1}, SubArray{T1}}, heights::Union{Vector{T2}, SubArray{T2}}, impdcs::Union{Vector{T2}, SubArray{T2}}, src_h::AbstractFloat, rec_h::AbstractFloat, nfreq::Int64, freqs::Vector{Int64} ) where {T1 <: Int64, T2 <: Float32}

    ihard = 0
    isoft = 0 
    flow_ress = [0.0, 0.0]
    atten = 0.0
    attenuations = Vector{Float64}()
    profile = Vector{Tuple{Float64, Float64}}()
    flowpr = Vector{Tuple{Float64, Float64}}()
    ignd = Vector{Int64}()


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
    points[3] = ( points[2][1] - 1 + res[1], res[2] )

    hillxz = [ profile[1] ]
    klocs = [1]
    nmm = 1
    for i in 1:3
        push!( hillxz, profile[ points[i][1] ] )
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
        dx, dy = hillxz[3] .- hillxz[1]
        hillxz[2] = hillxz[1] .+ 0.1 .* (dx, dy)
    end

    if hillxz[4][1] == hillxz[5][1]
        dx, dy = hillxz[5] .- hillxz[3]
        hillxz[4] = hillxz[3] .+ 0.1 .* (dx, dy)
    end

 # ================================================= Profile and Supporting Stuff Established ============================================================
 
  # ----------------------------------------------- Level Model ----------------------------------------------------- 
    
    if nmm <= 2
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
                duml, atten = egal( d/2, d/2, src_h, rec_h, flow_ress[2], flow_ress[2], 0.0, 0.0, 0.0, freqs[j], duml )
                atten = attenh
            end
            if isoft > 0
    # NON MI E' CHIARO IL SENSO DI "duml"
                duml, attens = egal( d/2, d/2, src_h, rec_h, flow_ress[1], flow_ress[1], 0.0, 0.0, 0.0, freqs[j], duml )
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

  # ------------------------------------------------- Hill Model ---------------------------------------------------- 
    
    zz2 = 0.0
    zcrit = 0.0
    if nmm == 3
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
        # Set up the source and receiver locations, normal to the corresponding plateaus
        cosθ = 1.0
        sinθ = 0.0
        Δx, Δz = hillxz[2] .- hillxz[1]
        if Δx > 0
            θ = atan(Δz, Δx)
            cosθ = cos(θ)
            sinθ = sin(θ)
        end

        # Source location is hs above the start of the terrain cut
        srcloc = [ hillxz[1] .+ (0, src_h) ]
        # Reflect the original source image
        push!( srcloc, srcloc[1] .+ (2*src_h*cosθ).*(sinθ, -cosθ) )

        cosθ = 1.0
        sinθ = 0.0
        Δx, Δz = hillxz[5] .- hillxz[4]
        if Δx > 0
            θ = atan(Δz, Δx)
            cosθ = cos(θ)
            sinθ = sin(θ)
        end

        # Right over the end of the cut receiver
        recloc = [ hillxz[5] .+ (0, rec_h) ]
        # reflect the original receiver image
        push!( recloc, recloc[1] .+ (2*rec_h*cosθ).*(sinθ, -cosθ) )

        for j in 1:nfreq
            if ihard > 0
                attenh = bakkernn( hillxz, srcloc, recloc, flow_ress[2], flow_ress[2], flow_ress[2], freqs[j] )
                atten = attenh
            end
            if isoft > 0
                attens = bakkernn( hillxz, srcloc, recloc, flow_ress[1], flow_ress[1], flow_ress[1], freqs[j] )
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

  # ------------------------------------------------ Valley Model --------------------------------------------------- 

    if nmm == 3 && zz2 < zcrit
        diff1 = profile[ klocs[2] ] .- profile[ klocs[1] ]
        diff2 = profile[ klocs[3] ] .- profile[ klocs[2] ]

        α1 = atan( diff1[2], diff1[1] )
        α2 = atan( diff2[2], diff2[1] )
        α = α2 - α1 + π
        # d0, d1 = @. √sum( [diff1^2, diff2^2] )
        d0 = @. √( diff1[1]^2 + diff1[2]^2 )
        d1 = @. √( diff2[1]^2 + diff2[2]^2 )

        for j in 1:nfreq
            attenh = 0.0
            attens = 0.0
            if ihard > 0
                attenh = dal( d0, d1, src_h, rec_h, α, flow_ress[2], flow_ress[2], freqs[j] )
                atten = attenh
            end
            if isoft > 0
                attens = dal( d0, d1, src_h, rec_h, α, flow_ress[1], flow_ress[1], freqs[j] )
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

    return attenuations
end

# Based on "https://www.geeksforgeeks.org/dda-line-generation-algorithm-computer-graphics/"
"""
Digital Differential Analyzer, for line rasterization
"""
function DDA( dtm::GeoArrays.GeoArray{Float32}, r0::Integer, c0::Integer, rn::Integer, cn::Integer )
    Δr = rn - r0
    Δc = cn - c0
    steps = max( abs(Δr), abs(Δc) )
    r_inc = Δr / steps
    c_inc = Δc / steps
    r = Float64(r0)
    c = Float64(c0)
    heigths_profile = Vector{Float32}()
    coords_profile = Vector{Tuple{Float64, Float64}}()
    for i in 1:steps
        rint, cint = round.(Int64, [r, c])
        push!( heigths_profile, dtm[rint, cint][1] )
        push!( coords_profile, Tuple{Float64, Float64}(ga.coords( map, [rint, cint])) )
        r += r_inc
        c += c_inc
    end
    return heigths_profile, coords_profile
end
#=  VERSIONE CON IL RASTER DELLE IMPEDENZE
function DDA( dtm::GeoArrays.GeoArray{Float32}, terrain_impedences::GeoArrays.GeoArray{Float32}, r0::Integer, c0::Integer, rn::Integer, cn::Integer )
    Δr = rn - r0
    Δc = cn - c0
    steps = max( abs(Δr), abs(Δc) )
    r_inc = Δr / steps
    c_inc = Δc / steps
    r = Float64(r0)
    c = Float64(c0)
    heigths_profile = Vector{Float32}()
    impedences_profile = Vector{Float32}()
    coords_profile = Vector{Tuple{Float64, Float64}}()
    for i in 1:steps
        rint, cint = round.(Int64, [r, c])
        push!( heigths_profile, dtm[rint, cint][1] )
        push!( impedences_profile, terrain_impedences[rint, cint][1] )
        push!( coords_profile, Tuple{Float64, Float64}(ga.coords( map, [rint, cint])) )
        r += r_inc
        c += c_inc
    end
    return heigths_profile, impedences_profiles, coords_profile
end
=#

function ground_loss( dtm::GeoArrays.GeoArray{Float32}, x0::Float64, y0::Float64, frequency::Int64, noData::Float32 )
    # Cell dimensions
    Δx::Float64, Δy::Float64 = ( ga.coords( dtm, size(dtm)[1:2] ) .- ga.coords( dtm, [1,1] ) ) ./ size(dtm)[1:2]
    # Source cell
    r0, c0 = ga.indices( dtm, [x0, y0] )
    # Source height
    h0::Float32 = dtm[r0, c0][1]
    # Maximum radius according to transmission loss
    max_radius = ceil(10^(dB / 20))
    # Number of cell fitting the radius of the area of effect of the sound
    cell_num = ceil( Int64, max_radius / max(Δx, Δy) )
    # Limits of the area
    row_begin = r0 - cell_num
    row_end = r0 + cell_num
    col_begin = c0 - cell_num
    col_end = c0 + cell_num
    # Trivial profiles
    # Profile containing the heights
    heights_profiles = push!(
        Vector{Vector{Float32}}(),
        dtm[ r0, c0:-1:col_begin ], # x axis first half
        dtm[ r0, c0:col_end ], # x axis second half
        dtm[ r0:-1:row_begin, c0 ], # y axis first half
        dtm[ r0:row_end, c0 ], # y axis second half
        [ dtm[r0+i, c0-i][1] for i in 0:-1:-cell_num ], # first diagonal first half
        [ dtm[r0+i, c0-i][1] for i in 0:cell_num ], # first diagonal second half
        [ dtm[r0+i, c0+i][1] for i in 0:-1:-cell_num ], # second diagonal first half
        [ dtm[r0+i, c0+i][1] for i in 0:cell_num ] # second diagonal second half
    )
    # Profile containing the positions
    coords_profiles = push!(
        Vector{Vector{Tuple{Float64, Float64}}}(),
        # x and y axes
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0, col])) for col in c0:-1:col_begin ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0, col])) for col in c0:col_end ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [row, c0])) for row in r0:-1:row_begin ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [row, c0])) for row in r0:row_end ],
        # diagonals
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0-i])) for i in 0:-1:-cell_num ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0-i])) for i in 0:cell_num ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0+i])) for i in 0:-1:-cell_num ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0+i])) for i in 0:cell_num ]
    )
    # Vector containing the indexes of the points on first quadrant of the border of the area of interest
    endpoints = vcat(
        [ (row_begin, col) for col in c0+1:col_end-1 ],
        [ (row, col_end) for row in row_begin+1:r0 ]
    )
    r, c = (0, 0)
    dists = nothing
    # Matrix with the resulting intenisty levels on the area of interest
    intenisty_matrix = fill( noData, row_end - row_begin + 1, col_end - col_begin + 1 )


 # PLACEHOLDER PER I PROFILI OTTENUTI DAL RASTER DELLE IMPEDENZE
    impdcs = zeros(Float32, length(heights_profiles[1]))


    # Compute and insert the values for the trivial profiles( x and y axes and diagonal profiles )
    @inbounds for (heights, coords) in zip(heights_profiles, coords_profiles)
        # Vector holding the distances from the source to the current point
        dists = [ distance((x0, y0), coords[1]) ]
        for j in 2:length(heights)
            # Add distance of the current point from source to the distances vector
            push!( dists, distance((x0, y0), coords[j]) )
            # Current cell
            r, c = ga.indices( dtm, [ coords[j]... ] )
            # If the cell should have a value
            if dtm[r, c] != noData
                intenisty_matrix[ ( (r, c) .- (row_begin, col_begin) .+ 1 )... ] = onCut( view(dists, 1:j), view(heights, 1:j), view(impdcs, 1:j), h0, heights[j], 1, frequency )[end]
            end
        end
    end
    # Compute the values for the remainder of the area
    @inbounds for point in endpoints
        for α in [0, 90, 180, 270]
            # Arrays of the heigths and the respective coordnates
            heights, impedences, coords = DDA( dtm, r0, c0, rotate_point( point..., r0, c0, α )... )
            # Compute array of the distances of each point of the profile from the source
            dists = [ distance((x0, y0), coords[1]) ]
            # Array of the resulting attenuations for each point of a single profile
            for j in 2:length(heights)
                push!( dists, distance((x0, y0), coords[j]) )
                r, c = ga.indices( dtm, [ coords[j]... ] )
                if dtm[r, c] != noData
                    intenisty_matrix[ ( (r, c) .- (row_begin, col_begin) .+ 1 )... ] = onCut( view(dists, 1:j), view(heights, 1:j), view(impdcs, 1:j), h0, heights[j], 1, frequency )[end]
                end
            end
        end
    end
    return intenisty_matrix 
end

# function noise_level( dtm::GeoArrays.GeoArray{Float32}, terrain_impedence::GeoArrays.GeoArray{Float32}, x0::Float64, y0::Float64, relative_humidity::Float64, temperature_k::FLoat64, frequencies::Union{Int64, Vector{Int64}}, noData::Float32 )
function noise_level( dtm::GeoArrays.GeoArray{Float32}, x0::Float64, y0::Float64, relative_humidity::Float64, temperature_k::FLoat64, frequency::Int64, noData::Float32 )
    # Cell dimensions
    Δx::Float64, Δy::Float64 = ( ga.coords( dtm, size(dtm)[1:2] ) .- ga.coords( dtm, [1,1] ) ) ./ size(dtm)[1:2]
    # Source cell
    r0, c0 = ga.indices( dtm, [x0, y0] )
    # Source height
    h0::Float32 = dtm[r0, c0][1]
    # Maximum radius according to transmission loss
    max_radius = ceil(10^(dB / 20))
    # Number of cell fitting the radius of the area of effect of the sound
    cell_num = ceil( Int64, max_radius / max(Δx, Δy) )
    # Limits of the area
    row_begin = r0 - cell_num
    row_end = r0 + cell_num
    col_begin = c0 - cell_num
    col_end = c0 + cell_num
    # Trivial profiles
    # Profile containing the heights
    heights_profiles = push!(
        Vector{Vector{Float32}}(),
        dtm[ r0, c0:-1:col_begin ], # x axis first half
        dtm[ r0, c0:col_end ], # x axis second half
        dtm[ r0:-1:row_begin, c0 ], # y axis first half
        dtm[ r0:row_end, c0 ], # y axis second half
        [ dtm[r0+i, c0-i][1] for i in 0:-1:-cell_num ], # first diagonal first half
        [ dtm[r0+i, c0-i][1] for i in 0:cell_num ], # first diagonal second half
        [ dtm[r0+i, c0+i][1] for i in 0:-1:-cell_num ], # second diagonal first half
        [ dtm[r0+i, c0+i][1] for i in 0:cell_num ] # second diagonal second half
    )
    # Profile containing the positions
    coords_profiles = push!(
        Vector{Vector{Tuple{Float64, Float64}}}(),
        # x and y axes
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0, col])) for col in c0:-1:col_begin ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0, col])) for col in c0:col_end ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [row, c0])) for row in r0:-1:row_begin ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [row, c0])) for row in r0:row_end ],
        # diagonals
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0-i])) for i in 0:-1:-cell_num ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0-i])) for i in 0:cell_num ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0+i])) for i in 0:-1:-cell_num ],
        [ Tuple{Float64, Float64}(ga.coords(dtm, [r0+i, c0+i])) for i in 0:cell_num ]
    )
    # Vector containing the indexes of the points on first quadrant of the border of the area of interest
    endpoints = vcat(
        [ (row_begin, col) for col in c0+1:col_end-1 ],
        [ (row, col_end) for row in row_begin+1:r0 ]
    )
    r, c = (0, 0)
    dists = nothing
    # Matrix with the resulting intenisty levels on the area of interest
    intenisty_matrix = fill( noData, row_end - row_begin + 1, col_end - col_begin + 1 )


 # PLACEHOLDER PER I PROFILI OTTENUTI DAL RASTER DELLE IMPEDENZE
    impdcs = zeros(Float32, length(heights_profiles[1]))


    # Compute and insert the values for the trivial profiles( x and y axes and diagonal profiles )
    @inbounds for (heights, coords) in zip(heights_profiles, coords_profiles)
        # Vector holding the distances from the source to the current point
        dists = [ distance((x0, y0), coords[1]) ]
        for j in 2:length(heights)
            # Add distance of the current point from source to the distances vector
            push!( dists, distance((x0, y0), coords[j]) )
            # Current cell
            r, c = ga.indices( dtm, [ coords[j]... ] )
            # If the cell should have a value
            if dtm[r, c] != noData
                intenisty_matrix[ ( (r, c) .- (row_begin, col_begin) .+ 1 )... ] =
                    # Loss of intensity due to propagation
                    transmission_loss( dists[j] ) -
                    # Loss of intensity due to the atmosphere
                    atmospheric_loss( dists[j], heights[j], relative_humidity, temperature_k, frequency ) -
                    # Loss of intenisty due to the terrain
                    onCut( view(dists, 1:j), view(heights, 1:j), view(impdcs, 1:j), h0, heights[j], 1, [dB] )[end]
            end
        end
    end
    # Compute the values for the remainder of the area
    @inbounds for point in endpoints
        for α in [0, 90, 180, 270]
            # Arrays of the heigths and the respective coordnates
            heights, impedences, coords = DDA( dtm, r0, c0, rotate_point( point..., r0, c0, α )... )
            # Compute array of the distances of each point of the profile from the source
            dists = [ distance((x0, y0), coords[1]) ]
            # Array of the resulting attenuations for each point of a single profile
            for j in 2:length(heights)
                push!( dists, distance((x0, y0), coords[j]) )
                r, c = ga.indices( dtm, [ coords[j]... ] )
                if dtm[r, c] != noData
                    intenisty_matrix[ ( (r, c) .- (row_begin, col_begin) .+ 1 )... ] =
                        # Loss of intensity due to propagation
                        transmission_loss( dists[j] ) -
                        # Loss of intensity due to the atmosphere
                        atmospheric_loss( dists[j], heights[j], relative_humidity, temperature_k, frequency ) -
                        # Loss of intenisty due to the terrain
                        onCut( view(dists, 1:j), view(heights, 1:j), view(impedences, 1:j), h0, heights[j], 1, [dB] )[end]
                end
            end
        end
    end
    return intenisty_matrix 
end






function do_noise( dtm::AbstractString, terrain_impedences::AbstractString, source, temperature_K::Float64, relative_humidity::Float64, frequencies )
    dtm_raster = replace( ga.read(dtm), missing => -9999.0f0 )
    terrain_raster = replace( ga.read(terrain_impedences), missing => -9999.0f0 )
    x0, y0 =

    data = noise_level( dtm_raster, terrain_raster, x0, y0, relative_humidity, temperature_K, frequencies )

    Functions.writeRaster()
end



















# ========================================================== TESTING =========================================================================================================

using BenchmarkTools

frequency = 110
dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
#   dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
dtm = replace( ga.read(dtm_file), missing => -9999.0f0 )
#   x0, y0 = (726454.9302346368, 5.025993899219433e6)
#       x0, y0 = (11.930065824163105,45.425861311724816)
x0, y0 = ga.coords(dtm, [4500, 5700])

GC.gc()

@code_warntype ground_loss(dtm, x0, y0, frequency)

@time ground_loss(dtm, x0, y0, frequency)

@btime ground_loss(dtm, x0, y0, frequency)

mat = noise_level(dtm, x0, y0, 0.2, 293.15, frequency)

cmat = map( x -> x = x == -9999.0f0 ? -30.0f0 : x, mat )
cols, rows = size(cmat)
# Risultato
heatmap( 1:cols, 1:rows, cmat )
# DTM della stessa area del risultato
heatmap( 1:637, 1:637, dtm.A[ 4500-(637÷2):4500+(637÷2), 5700-(637÷2):5700+(637÷2), 1 ] )


#=
attenuations_points = []
for point in endpoints
    if point[2] == c0 || abs(point[1] - r0) == abs(point[2] - c0)
        continue
    end
    for α in [0, 90, 180, 270]
        rm, cm = rotate_point( point[1], point[2], r0, c0, α )
        # Arrays of the heigths and the respective coordnates
        heights, coords = DDA( dtm, r0, c0, rm, cm )
        # Compute array of the distances of each point of the profile from the source
        dists = map( p -> √( ( p[1] - x0 )^2 + ( p[2] - y0 )^2 ), coords )
        # Array of the resulting attenuations for each point of a single profile
        atten_point = []
        for j in 2:length(heights)
            atten = onCut(dists[1:j], heights[1:j], zeros(length(heights[1:j])), h0, heights[j], 1, [dB] )
            x, y = ga.indices( dtm, [ coords[j]... ] ) .- [ row_begin, col_begin ] 
            push!( atten_point, ( [x, y], [coords[j]...], atten ) )
        end
        push!( attenuations_points, atten_point )
    end
end

mat = Array{Any}( missing, row_end - row_begin, col_end - col_begin )
for ( profile, results ) in zip( coords_profiles, attenuations )
    for (coords, atten) in zip(profile, results)
        row, col = Int64.( ga.indices(dtm, [coords...]) .- [row_begin, col_begin] )
        if ismissing(mat[row, col]) || mat[row, col][3]  < atten[1]
            mat[row, col] = ( coords[1], coords[2], atten[1] )
        end
    end
end

# DA SISTEMARE
mat = Array{Any}( missing, row_end - row_begin, col_end - col_begin )
for points in attenuations_points
    for point in points
        if ismissing(mat[point[1]...])
            mat[point[1]...] = ( point[1][1], point[1][2], point[3][1] )
        else
            if mat[point[1]...][3] != point[3][1]
                #   println("Attennuazione diversa per lo stesso punto")
                #   println("$(mat[point[1]...]) e $point\n")
                if point[3][1] > mat[point[1]...][3]
                    mat[point[1]...] = ( point[1][1], point[1][2], point[3][1] )
                end
            end
        end
    end
end

plot( mat[1,2], seriestype=:scatter )
for point in mat
    if !ismissing(point)
        plot!(point, seriestype=:scatter)
    end
end
current()
=#


noData = -9999.f0
mat = fill( noData, row_end - row_begin, col_end - col_begin )
for point in endpoints
    for α in [0, 90, 180, 270]
        rm, cm = rotate_point( point[1], point[2], r0, c0, α )
        # Arrays of the heigths and the respective coordnates
        heights, coords = DDA( dtm, r0, c0, rm, cm )
        # Compute array of the distances of each point of the profile from the source
        dists = map( p -> √( ( p[1] - x0 )^2 + ( p[2] - y0 )^2 ), coords )
        # Array of the resulting attenuations for each point of a single profile
        for j in 2:length(heights)
            atten = onCut(dists[1:j], heights[1:j], zeros(length(heights[1:j])), h0, heights[j], 1, [dB] )[end]
            r, c = ga.indices( dtm, [ coords[j]... ] ) .- [ row_begin, col_begin ]
            if mat[r, c] == noData || atten > mat[r, c]
                mat[r, c] = atten
            end
        end
    end
end

cmat = map( x -> x = x == noData ? 0.0f0 : x, mat )
cols, rows = size(cmat)
heatmap( 1:cols, 1:rows, cmat, c=cgrad([:blue, :white, :yellow, :red]) )



cols, rows = size(dtm)
heatmap( 1:cols, 1:rows, mat, c=cgrad([:blue, :white, :yellow, :red]) )


#=
    #   x0, y0 = (710705.0, 5065493.0)
    x0, y0 = (726454.9302346368, 5.025993899219433e6)
    #   x0, y0 = (11.930065824163105,45.425861311724816)
    dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff"
    #   dtm_file = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
    dtm = Raster(dtm_file)
    # Coordinate dei punti
    x1, y1 = Functions.toCoords(dtm, 1, 1)
    xn, yn = Functions.toCoords( dtm, size(dtm)[1:2]... )
    # Dimensioni in metri di una cella
    Δx, Δy = (xn-x1, yn-y1) ./ size(dtm)[1:2]
    dB = 110 - 32
    h0 = dtm[x0, y0]
    max_radius = ceil(10^(dB/20))
    cell_num = ceil( Int64, max_radius/Δx )
    row_begin = x0 - max_radius
    row_end = x0 + max_radius
    col_begin = y0 + max_radius
    col_end = y0 - max_radius

    # Trivial profiles
    # Profile containing the heights
    heights_profiles = [
        dtm[ x0, y0:-Δy:col_begin ], # y axis first half
        dtm[ x0, y0:Δy:col_end ], # y axis second half
        dtm[ x0:-Δx:row_begin, y0 ], # x axis first half
        dtm[ x0:Δx:row_end, c0 ], # x axis second half
        [ dtm[x0 + (i * Δx), y0 - (i * Δy)] for i in 0:-1:-cell_num ], # first diagonal first half
        [ dtm[x0 + (i * Δx), y0 - (i * Δy)] for i in 0:cell_num ], # first diagonal second half
        [ dtm[x0 + (i * Δx), y0 + (i * Δy)] for i in 0:-1:-cell_num ], # second diagonal first half
        [ dtm[x0 + (i * Δx), y0 + (i * Δy)] for i in 0:cell_num ]  # second diagonal second half
    ]
    # Profile containing the positions
    coords_profiles = [
        # x and y axes
        [ Functions.toCoords(dtm, x0, col) for col in c0:-1:col_begin ],
        [ Functions.toCoords(dtm, x0, col) for col in c0:col_end ],
        [ Functions.toCoords(dtm, row, y0) for row in r0:-1:row_begin ],
        [ Functions.toCoords(dtm, row, y0) for row in r0:row_end ],
        # diagonals
        [ Functions.toCoords(dtm, x0 + (i * Δx), y0 - (i * Δy)) for i in 0:-1:-cell_num ],
        [ Functions.toCoords(dtm, x0 + (i * Δx), y0 - (i * Δy)) for i in 0:cell_num ],
        [ Functions.toCoords(dtm, x0 + (i * Δx), y0 + (i * Δy)) for i in 0:-1:-cell_num ],
        [ Functions.toCoords(dtm, x0 + (i * Δx), y0 + (i * Δy)) for i in 0:cell_num ]
    ]

    top_idxs = [ (row_begin, col) for col in  y0+Δy:Δy:col_end-Δy ]
    right_idxs = [ (row, col_end) for row in  row_begin+Δx:Δx:x0 ]
    # Vector containing the indexes of the points on first quadrant of the border of the area of interest
    endpoints = vcat( top_idxs, right_idxs )



    #=
    attenuations_points = []
    for point in endpoints
        if point[2] == c0 || abs(point[1] - r0) == abs(point[2] - c0)
            continue
        end
        for α in [0, 90, 180, 270]
            rm, cm = rotate_point( point[1], point[2], r0, c0, α )
            # Arrays of the heigths and the respective coordnates
            heights, coords = DDA( dtm, r0, c0, rm, cm )
            # Compute array of the distances of each point of the profile from the source
            dists = map( p -> √( ( p[1] - x0 )^2 + ( p[2] - y0 )^2 ), coords )
            # Array of the resulting attenuations for each point of a single profile
            atten_point = []
            for j in 2:length(heights)
                atten = onCut(dists[1:j], heights[1:j], zeros(length(heights[1:j])), h0, heights[j], 1, [dB] )
                x, y = ga.indices( dtm, [ coords[j]... ] ) .- [ row_begin, col_begin ] 
                push!( atten_point, ( [x, y], [coords[j]...], atten ) )
            end
            push!( attenuations_points, atten_point )
        end
    end

    mat = Array{Any}( missing, row_end - row_begin, col_end - col_begin )
    for ( profile, results ) in zip( coords_profiles, attenuations )
        for (coords, atten) in zip(profile, results)
            row, col = Int64.( ga.indices(dtm, [coords...]) .- [row_begin, col_begin] )
            if ismissing(mat[row, col]) || mat[row, col][3]  < atten[1]
                mat[row, col] = ( coords[1], coords[2], atten[1] )
            end
        end
    end


    # DA SISTEMARE
    mat = Array{Any}( missing, row_end - row_begin, col_end - col_begin )
    for points in attenuations_points
        for point in points
            if ismissing(mat[point[1]...])
                mat[point[1]...] = ( point[1][1], point[1][2], point[3][1] )
            else
                if mat[point[1]...][3] != point[3][1]
                    #   println("Attennuazione diversa per lo stesso punto")
                    #   println("$(mat[point[1]...]) e $point\n")
                    if point[3][1] > mat[point[1]...][3]
                        mat[point[1]...] = ( point[1][1], point[1][2], point[3][1] )
                    end
                end
            end
        end
    end

    plot( mat[1,2], seriestype=:scatter )
    for point in mat
        if !ismissing(point)
            plot!(point, seriestype=:scatter)
        end
    end
    current()
    =#



    xm, ym = endpoints[273]
    res_h, res_c = DDA(dtm, x0, y0, xm, ym, Δx, Δy)

    plot([(x0, y0), (xm, ym)])
    plot!(res_c)

    noData = -9999.f0
    mat = fill( noData, round(Int64, (row_end - row_begin) / Δx), round(Int64, (col_end - col_begin) / Δy) )
    for point in endpoints
        for α in [0, 90, 180, 270]
            rm, cm = rotate_point( point[1], point[2], r0, c0, α )
            # Arrays of the heigths and the respective coordnates
            heights, coords = DDA( dtm, r0, c0, rm, cm )
            # Compute array of the distances of each point of the profile from the source
            dists = map( p -> √( ( p[1] - x0 )^2 + ( p[2] - y0 )^2 ), coords )
            # Array of the resulting attenuations for each point of a single profile
            for j in 2:length(heights)
                atten = onCut(dists[1:j], heights[1:j], zeros(length(heights[1:j])), h0, heights[j], 1, [dB] )[end]
                r, c = ga.indices( dtm, [ coords[j]... ] ) .- [ row_begin, col_begin ]
                if mat[r, c] == noData || atten > mat[r, c]
                    mat[r, c] = atten
                end
            end
        end
    end

    #   cmat = deepcopy(mat)
    cmat = map( x -> x = x == noData ? 0.0 : x, mat )
    cols, rows = size(cmat)
    heatmap( 1:cols, 1:rows, cmat, c=cgrad([:blue, :white, :yellow, :red]) )



    cols, rows = size(dtm)
    heatmap( 1:cols, 1:rows, mat, c=cgrad([:blue, :white, :yellow, :red]) )
=#

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

    size = round(Int64, √length(centers))
    start = size * floor(Int64, size/2) + 1
    x_axis = centers[ start:start+size-1 ]
    plot!(x_axis, c=:blue )

    start = ceil(Int64, size/2)
    y_axis = centers[[ start + size*i for i in 0:size-1 ]]
    plot!(y_axis, c=:red )

    diagonal1 = centers[[ size + (size-1)i for i in 0:size-1 ]]
    plot!(diagonal1, c=:yellow)

    diagonal2 = centers[[ 1 + (size+1)i for i in 0:size-1 ]]
    plot!(diagonal2, c=:purple)

    is = collect(x0:size)
    x0 = y0 = Int64(round(size/2))
    r = floor(Int64, size/2)
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