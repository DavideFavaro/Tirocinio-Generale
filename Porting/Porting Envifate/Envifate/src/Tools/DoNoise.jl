module DoNoise
"""
"""


using GeoStats


Base.:-( x::Tuple, y::Tuple ) = ( x[1] - y[1], x[2] - y[2] )


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
    #   frO = ( p_atm_pressure * ( (24 + 4.04e04) * humidity ) * (0.02 + humidity) ) / (0.391 + humidity)
    frO = p_atm_pressure * (
                               24 + ( 
                                        4.04e04 * humidity * (
                                                                 (0.02 + humidity) / (0.391 + humidity)
                                                             )
                                    )
                           )

    # Calculate relaxation frequency of N (equation 4)
    frN = p_atm_pressure * √t_tr * (
                                       9 + (
                                               280 * humidity * ℯ^( -4.170 * (
                                                                                 t_tr^(-1/3) - 1
                                                                             )
                                                                  )
                                           )
                                   )

    # Calculate alpha (equation 5)
    term1 = 1.84*10^(-11) * p_atm_pressure^(-1) * √t_tr
    #   term2 = t_tr^(-2.5) * ( 0.01275 * ℯ^(-2239.1 / temp_k) * ( frO / (frO^2 + freq^2) ) )
    term2 = t_tr^(-2.5) * (
                              0.01275 * ℯ^(-2239.1 / temp_k) / (
                                                                   frO + (freq^2 / frO)
                                                               )
                          )
    #   term3 = 0.1068 * ℯ^(-3352 / temp_k) * ( frN / (frN^2 + freq^2) )
    term3 = t_tr^(-2.5) * (
                              0.1068 * ℯ^(-3352 / temp_k) / (
                                                                frN + (freq^2 / frN)
                                                            )
                          ) 
    #   α = 8.686 * (frequency^2) * ( term1 + term2 + term3 )
    α = frequency^2 * (term1 + term2 + term3)

    return α * r / 100
end






function minmax( profile::Vector{Tuple{Float64, Float64}}, rel_h_src::Real=0.0, rel_h_rec::Real=0.0 )
    # Define parameters A and B such that height of the source-receiver line
    # is given by z = Ax + B.
    x1 = profile[1][1]
    z1 = profile[1][2] + rel_h_src
    xn = profile[end][1]
    zn = profile[end][2] + rel_h_rec
    a = (zn - z1) / (xn - x1)
    b = z1 - a*x1

    # Look for the max and min
    prf = map( row -> ( row[1], row[2] - ( a*row[1] + b ) ), profile )
    sort!( prf, by=(row) -> row[2] )
    return ( prf[1], prf[end] ) 
end

function delbaz( freq::Real, flow_res::Real )::Complex
    dumr = 1.0 + 9.08 * (flow_res/freq)^0.75
    dumi = -11.9 * (flow_res/freq)^0.73
    return complex(dumr, dumi)
end

function subw2( a::Real, b::Real, c::Real, d::Real, w::Complex )::Complex
    an = (a^2 - b^2  - d)^2 + (2 * a * b)^2
    ir = c * (a^2 + b^2 + d)
    wr = b * ir / an + w.re
    wi = a * ir / an + w.im
    return complex( wr, wi )
end

function ww2(t::Complex)::Complex
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
            w = subw2(a, b, c, d, w)    
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

function qq2( d::Real, src_h::Real, rec_h::Real, freq::Real, z::Complex )::Complex
    k = 2.0 * π * freq/340.0
    #   r1 = √( d^2 + (src_h - rec_h)^2 )
    r = √( d^2 + (src_h + rec_h)^2 )
    c = (src_h + rec_h) / r
    
    n = (z.re * c + 1.0)^2 + (z.im * c)^2
    rr = ( ( c * abs(z) )^2 - 1 ) / n
    ri = 2.0 * z.im * c / n

    nyr = z.re / abs2(z)
    nyi = z.im / abs2(z)

    dumr = √(k * r / 4.0) * (nyr + c + nyi)
    dumi = √(k * r / 4.0) * (-nyr - c + nyi)
    t = complex(dumr, dumi)

    w = ww2(-t)

    fr = 1.0 + √π * real(t*w)
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
    arg8  = transition_height ? o turbulence? o cosa?
    arg9  = ? 
    arg10 = ?
    arg11 = ?
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

function oncut( dists::AbstractVector, heights::AbstractVector, impdcs::AbstractVector, src_h::Real, rec_h::Real, nfreq::Int64, freqs::AbstractVector )
 #=
    hgts(i) ==> points[i][1]
    locs(i) ==> points[i][2]
 =#
    ihard = isoft = 0 
    flohard = flosoft = 0.0
    profile = flowpr = ignd = []

    for i in 1:length(dists)
        push!( profile, (dists[i], heights[i]) )
        push!( flowpr, (dists[i], impdcs[i]) )
        if flowpr[i][2] <= 1000.0
            push!( ignd, 0 )
            isoft += 1
            flosoft += flowpr[i][2]
        else
            push!( ignd, 1 )
            ihard += 1
            flohard += flowpr[i][2]
        end
    end
    flosoft /= max(isoft)
    flohard /= max(ihard)
    
    flosoft == 0.0 && flosoft = 200.0
    floHard == 0.0 && floHard = 10^6


    points = [ (0.0, 0.0), (0.0, 0.0), (0.0, 0.0) ]
    points[2] = minmax( profile, rel_h_src, rel_h_rec )[2]
    points[1] = minmax( profile[ 1:points[2][1] ], rel_h_src, 0.0 )[1]
    ntemp = length(profile) - points[2][1] + 1
    #   call maxmin(prof(1,locs(2)),ntemp,0.,hrec,hgtmm,locmm)
    # Non so cosa voglia dire `prof(1,locs(2))`
    res = minmax( profile[ point[2][1] ][1], 0.0, rel_h_rec  )
    points[3][1] = points[2][1] + res[1] - 1
    points[3][2] = res[2] 


    hillxz[1] = [ (profile[1][1], profile[1][2]) ]
    klocs = [1]
    nmm = 1
    for i in 1:3
        hill = (profile[ points[i][2] ][1], profile[ points[i][2] ][2])
        push!( hillxz, hill )
        if point[i][2] != klocs[nmm]
            nmm += 1
            push!( klocs, points[i][2] )
        end
    end

    push!( hillxz, (profile[end][1], profile[end][2]) )
    if klocs[nmm] < length(profile)
        nmm += 1
        push!( klocs, length(profile) )
    end


    if points[1][2] == points[2][2]
        hillxz[2] = hillxz[1]
    end

    if points[3][2] == points[2][2]
        # hillxz[4] = hillxz[5]
        hillxz[end-1] = hillxz[end]
    end

    if hillxz[1] == hillxz[2]
        dx, dy = hillxz[3] - hillxz[1]
        hillx[2] = ( hillxz[1][1] + 0.1*dx, hillxz[1][2] + 0.1*dy  )
    end

    if hillxz[1] == hillxz[2]
        dx, dy = hillxz[5] - hillxz[3]
        hillx[4] = ( hillxz[3][1] + 0.1*dx, hillxz[3][2] + 0.1*dy  )
    end


    if nmm <= 2
        ax = profile[1][1]
        ox = profile[end][1]

        #   dist = √( (ax-ox)^2 + (ay-oy)^2 )
        # `ay` e `oy` sono uguali a zero quindi la funzione diventa `√( (ax-ox)^2 )`
        d = abs(ax - ox)

        for j in 1:nfreq
            dumf = freq[j]
            if ishard > 0
                duml, attenh = egal( d/2, d/2, src_h, rec_h, flohard, flohard, vl, 0.0, 0.0, dumf, duml )
                atten = attenh
            end
            if isoft > 0
                duml, attens = egal( d/2, d/2, src_h, rec_h, flosoft, flosoft, vl, 0.0, 0.0, dumf, duml )
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( dists, ignd, src_h, rec_h, attens, attenh )
            end
        end
    elseif nmm == 3
        dist = profile[klocs[3]][1] - profile[klocs[1]][1]
        dsl = profile[klocs[2]][1] - profile[klocs[1]][1]
        zz1 = profile[klocs[1]][2]
        zz2 = profile[klocs[2]][2]
        zz3 = profile[klocs[3]][2]
        zcrit = zz1 + (zz3-zz1) * dsl / dist
    elseif nmm > 3 || ( nmm == 3 && zz2 >= zcrit ) 
        costh = 1
        sinth = 0
        delx = hillxz[2][1] - hillxz[1][1]
        if delx > 0
 # NON SO QUALE SIA IL CORRISPETTIVO DI ATAN2()
            thet = atan( )
            costh = cos(thet)
            sinth = sin(thet)
    else
    end 
end









function create_Grid( xs::Real, ys::Real, intensity::Real )
    # maximum ideal radius for the diffusion of a sound with given intensity  
    max_radius = 10^(intensity/20)
    # Coordinates of the first point (upper left) of the grid containing the radius 
    x0 = xs - max_radius 
    y0 = ys + max_radius

    grid = CartesianGrid( (x0, ys-max_radius), (xs+max_radius, y0), (200.,200.) )

    return grid
end



end # module