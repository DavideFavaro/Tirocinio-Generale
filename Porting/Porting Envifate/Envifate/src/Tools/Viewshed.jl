module Viewshed


# Rotate a point P around a center C by angle θ
rotate_x( xp, yp, xc, yc, θ ) = round( Int64, (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) + xc )
rotate_y( xp, yp, xc, yc, θ ) = round( Int64, (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) + yc )
rotate_point( xp, yp, xc, yc, θ ) = θ == 0 ? (xp, yp) : ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )

# Given a raster, create the profiles by taking all the cells between the one at `(c0, r0)` and the one at `(rn, cn)`  
function DDA( dtm, r0::Number, c0::Number, rn::Number, cn::Number )
    Δr = rn - r0
    Δc = cn - c0
    steps = max( abs(Δr), abs(Δc) )
    r_inc = Δr / steps
    c_inc = Δc / steps
    r_end = r0 + r_inc * steps
    c_end = c0 + c_inc * steps 
    heigths_profile = []
    coords_profile = []
    for (r, c) in zip( r0:r_inc:r_end, c0:c_inc:c_end ) 
        rint, cint = round.(Int64, [r, c])
        push!( heigths_profile, dtm[rint, cint][1] )
        push!( coords_profile, Tuple(GeoArracs.coords( dtm, [rint, cint])) )
    end
    return heigths_profile, coords_profile
end


function generate_profile( x0, y0, xn, yn )
    r0, c0 = toIndexes(dtm, x0, y0)
    rn, cn = toIndexes(dtm, xn, yn)
    return DDA(dtm, r0, c0, rn, cn)
end

# Given the number of points in a profile, the indexes of the source and the raster, generate 360 profiles 
function generate_profiles( points_number, x0, y0, dtm )
    values_profiles = []
    indexes_profiles = []
    r0, c0 = toIndexes(dtm, x0, y0)
    rn = r0+points_number
    for α in [0, 90, 180, 270]
        for β in 0:89
            rn, cn = rotate_point(rn, c0, r0, c0, α+β)
            values, indexes = DDA(dtm, r0, c0, rn, cn)
            push!(values_profiles, values)
            push!(indexes_profiles, indexes)
        end
    end
    return values_profiles, indexes_profiles
end

# Given a profile, return only the visible points from the source, if `result == :visible`, all the invisible ones otherwise  
function veiwshed( profile::AbstractVector, result::Symbol=:visible )::AbstractVector
    if result ∉ [:visible, :invisible]
        throw(DomainError(result, "result must either be :visible or :invisible"))
    end
    if result == :visible
        visible = [ profile[1], profile[2] ]
        slope = ( profile[2][2] - profile[1][2] ) / ( profile[2][1] - profile[1][1] )
        for i in 3:length(profile)
            new_slope = ( profile[i][2] - profile[1][2] ) / ( profile[i][1] - profile[1][1] )
            if new_slope >= slope
                push!( visible, profile[i] )
                slope = new_slope
            end
        end
        return visible
    else
        invisible = []
        for i in 3:length(profile)
            new_slope = ( profile[i][2] - profile[1][2] ) / ( profile[i][1] - profile[1][1] )
            if new_slope < slope
                push!( nonvisible, profile[i] )
                slope = new_slope
            end
        end
        return invisible
    end
end