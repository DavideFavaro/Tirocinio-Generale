module Viewshed
"""
Module with functions for the computation of profiles and viewshed
"""


#=
    FORSE SI PUO' MODIFICARE NELL'OTTICA DI USARE Rasters.jl E IN GENERALE SI POTREBBE METTERE UN PO' A POSTO
=#


"""
    rotate_x( xp::Number, yp::Number, xc::Number, yc::Number, θ::Number )::Int64

Find the cartesian `x` coordinate of a point in `xp` after a rotation of `θ` degrees around a point in `xc`.

Used as part of the `rotate_point` function.
"""
rotate_x( xp::Number, yp::Number, xc::Number, yc::Number, θ::Number )::Int64 = round( Int64, (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) + xc )

"""
    rotate_x( xp::Number, yp::Number, xc::Number, yc::Number, θ::Number )::Int64

Find the cartesian `y` coordinate of a point in `yp` after a rotation of `θ` degrees around a point in `yc`.

Used as part of the `rotate_point` function.
"""
rotate_y( xp::Number, yp::Number, xc::Number, yc::Number, θ::Number ) = round( Int64, (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) + yc )

"""
    rotate_point( xp::Number, yp::Number, xc::Number, yc::Number, θ::Number )::Tuple{Int64, Int64}

Rotate a point of cartesian coordinates (`xp`, `yp`) around a point (`xc`, `yc`) by `θ` degrees.
"""
rotate_point( xp::Number, yp::Number, xc::Number, yc::Number, θ::Number ) = θ == 0 ? (xp, yp) : ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )



# SAREBBE PIU' CORRETTO METTERE dtm::Matrix{T} E DENTRO LA TUPLE DI RITORNO COME PRIMO ELEMENTO Vector{T}
"""
    getprofile( dtm::AbstractArray, r0::Integer, c0::Integer, rn::Integer, cn::Integer )::Tuple{ Vector{T}, Vector{ Tuple{Real, Real} } } where {T <: Number}

Return a `Tuple` of `Vectors` containing the values and the coordinates of all cells on the line between the point at indexes (`r0`, `c0`) and the one at (`rn`, `cn`) of matrix `dtm`.

The function is an implementation of the DDA (Digital Differential Analyser) algorithm for line rasterization.
""" 
function getprofile( dtm::AbstractArray, r0::Integer, c0::Integer, rn::Integer, cn::Integer )
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
        push!( coords_profile, Tuple(GeoArrays.coords( dtm, [rint, cint])) )
    end
    return heigths_profile, coords_profile
end

"""
    generate_profile( dtm::AbstractArray, x0::Real, y0::Real, xn::Real, yn::Real )::Tuple{ Vector{T}, Vector{ Tuple{Real, Real} } } where {T <: Number}

Return the profile, formed by a `Vector` of heights and one of coordinates, that spans from the point of coordinates (`x0`, `y0`) to the one at (`xn`, `yn`).
"""
function generate_profile( dtm::AbstractArray, x0::Real, y0::Real, xn::Real, yn::Real )
    r0, c0 = toIndexes(dtm, x0, y0)
    rn, cn = toIndexes(dtm, xn, yn)
    return getprofile(dtm, r0, c0, rn, cn)
end

# Given the number of points in a profile, the indexes of the source and the raster, generate 360 profiles
"""
    generate_profiles( dtm, x0::Real, y0::Real, points_number::Integer )::Tuple{ Vector{ Vector{T} }, Vector{ Vector{ Tuple{Integer, Integer} } } } where {T <: Number}

Given a raster `dtm`, a source point of coordinates (`x0`, `y0`) and `points_number` number of desired points for each profile,
return two `Vector`s containing, respectively, the `Vector`s of the heights and the corresponding indexes of the 360 profiles of `dtm`,
obtained as radii with the point as center.
"""
function generate_profiles( dtm::AbstractArray, x0::Real, y0::Real, points_number::Integer )
    values_profiles = []
    indexes_profiles = []
    r0, c0 = toIndexes(dtm, x0, y0)
    rn = r0+points_number
    for α in [0, 90, 180, 270]
        for β in 0:89
            rn, cn = rotate_point(rn, c0, r0, c0, α+β)
            values, indexes = getprofile(dtm, r0, c0, rn, cn)
            push!(values_profiles, values)
            push!(indexes_profiles, indexes)
        end
    end
    return values_profiles, indexes_profiles
end



# Given a profile, return only the visible points from the source, if `result == :visible`, all the invisible ones otherwise  
"""
    veiwshed( profile::AbstractVector, result::Symbol=:visible )::AbstractVector

Return a `Vector` of all the visible points of `profile`, if `result` is set to `:visible`, or of all the hidden points if set to `:hidden`.

The functions assumes the origin point, from which the viewshed is computed, to be in the first cell of the array.
"""
function veiwshed( profile::AbstractVector, result::Symbol=:visible )::AbstractVector
    if result ∉ [:visible, :hidden]
        throw(DomainError(result, "result must either be `:visible` or `:hidden`"))
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
        hidden = []
        for i in 3:length(profile)
            new_slope = ( profile[i][2] - profile[1][2] ) / ( profile[i][1] - profile[1][1] )
            if new_slope < slope
                push!( nonvisible, profile[i] )
                slope = new_slope
            end
        end
        return hidden
    end
end