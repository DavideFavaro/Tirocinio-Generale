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
    function getprofile( dtm::AbstractArray, r0::Integer, c0::Integer, rn::Integer, cn::Integer, h0::Real )

Return a `Tuple` of `Vectors` containing the values and the coordinates of all cells on the line between the point at indexes (`r0`, `c0`) and the one at (`rn`, `cn`) of matrix `dtm`.
The algorithm will stop at the first cell with value greater than the sum of `h0` and the value of cell (`r0`, `c0`).

The function is an implementation of the DDA (Digital Differential Analyser) algorithm for line rasterization.
"""
function getprofile( dtm::AbstractArray, r0::Integer, c0::Integer, rn::Integer, cn::Integer, h0::Real )
    Δr = rn - r0
    Δc = cn - c0
    # Height of the light source
    th0 = dtm[r0, c0] + h0
    # Number of steps to reach the final point
    steps = max( abs(Δr), abs(Δc) )
    r_inc = Δr / steps
    c_inc = Δc / steps
    #   r_end = r0 + r_inc * steps
    #   c_end = c0 + c_inc * steps 
    heigths_profile = []
    coords_profile = []
    for (r, c) in zip( r0:r_inc:rn, c0:c_inc:cn ) 
        rint, cint = round.(Int64, [r, c])
        push!( heigths_profile, dtm[rint, cint][1] )
        push!( coords_profile, Tuple(GeoArrays.coords( dtm, [rint, cint])) )
        if dtm[rint, cint] >= th0
            break
        end
    end
    return heigths_profile, coords_profile
end

#= CREDO SI POSSANO TOGLIERE
    """
        generate_profile( dtm::AbstractArray, x0::Real, y0::Real, xn::Real, yn::Real )::Tuple{ Vector{T}, Vector{ Tuple{Real, Real} } } where {T <: Number}

    Return the profile, formed by a `Vector` of heights and one of coordinates, that spans from the point of coordinates (`x0`, `y0`) to the one at (`xn`, `yn`).
    """
    function generate_profile( dtm::AbstractArray, x0::Real, y0::Real, xn::Real, yn::Real )
        r0, c0 = toIndexes(dtm, x0, y0)
        rn, cn = toIndexes(dtm, xn, yn)
        return getprofile(dtm, r0, c0, rn, cn)
    end

    function generate_profile( dtm::AbstractArray, x0::Real, y0::Real, h0::Real )
        r0, c0 = toIndexes(dtm, x0, y0)
        return getprofile(dtm, r0, c0, h0)
    end
=#

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
    rn = r0 + points_number
    for α in [0, 90, 180, 270], β in 0:89
        rn, cn = rotate_point(rn, c0, r0, c0, α+β)
        values, indexes = getprofile(dtm, r0, c0, rn, cn)
        push!(values_profiles, values)
        push!(indexes_profiles, indexes)
    end
    return values_profiles, indexes_profiles
end

"""
Generate the profiles from the source point to the first point higher than the source (the one that will hide all the following) 
"""
function generate_profiles( dtm::AbstractArray, x0::Real, y0::Real, h0::Real )
    values_profiles = []
    indexes_profiles = []
    r0, c0 = toIndexes(dtm, x0, y0)
    rn = size(dtm, 1)
    for α in [0, 90, 180, 270], β in 0:89
        rn, cn = rotate_point(rn, c0, r0, c0, α+β)
        values, indexes = getprofile(dtm, r0, c0, rn, cn, h0)
        push!(values_profiles, values)
        push!(indexes_profiles, indexes)
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
        throw(DomainError(result, "result must either be `:visible` or `:hidden`."))
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






distance( x0, y0, x1, y1 ) = √( ( x1 - x0 )^2 + ( y1 - y0 )^2 )
distance( p0, p1 ) = √( ( p1[1] - p0[1] )^2 + ( p1[2] - p0[2] )^2 )

"""
    viewshed( dtm, x0::Real, y0::Real, h0::Real )

Return the visibility matrix containing boolean values indicating wether a cell is reached by the light generated by a source positioned at (`x0`, `y0`) and with relative height of `h0`.
The matrix will have the same size of `dtm`. 
"""
function viewshed( dtm, x0::Real, y0::Real, h0::Real )
    # Source cell
    r0, c0 = toIndexes(dtm, x0, y0)
    # Final point of the right horizontal profile 
    rm = size(dtm, 1)
    # Total height of the source accounting for terrain
    th0 = dtm[r0, c0] + h0
    # Visibility matrix
    vizmat = falses(size(dtm)...)
    # The source is always visible
    vizmat[r0, c0] = true
    # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent.
    for α in 0:89, β in [0, 90, 180, 270]
        # Final point of the profile of the line with agle `α + β`
        rn, cn = rotate_point( rm, c0, r0, c0, α + β )
        Δr, Δc = (rn - r0, cn - c0)
        # Values to add to row and column of a cell on the line to pass to another cell of the line
        r_inc, c_inc = (rn - r0, cn - c0) ./ max( abs(Δr), abs(Δc) )
        # Indexes of the first cell after the source cell
        r1, c1 = round.(Int64, (r0, c0) .+ c_inc)
        # The first cell after the source is always visible
        vizmat[r1, c1] = true
        # Slope between the source and the first cell
        slope = ( dtm[r1, c1] - dtm[r0, c0] ) / distance( toCoords.(Ref(dtm), r1, c1), (0, 0) )
        # Iterate over each cell of the profile on the line
        for (r, c) in zip(r1:r_inc:rn, c1:c_inc:cn)
            # Calculate the precise indexes of the cell
            rint, cint = round.(Int64, [r, c])
            # Skip cell already known to be visible
            !vizmat[rint, cint] && continue
            # Compute the slope of the new cell 
            new_slope = (dtm[rint, cint] - dtm[r0, c0]) / distance( toCoords.(Ref(dtm), r1, c1), (0, 0) )
            # If the new slope is greater than the original one the cell is visible
            if new_slope >= slope
                vizmat[rint, cint] = true
                slope = new_slope
            end
            # If the cell is higher than the source al cell beyond the current ne will be hidden
            dtm[rint, cint] > th0 && break
        end
    end
    return vizmat
end

#= References
    function getprofile( dtm::AbstractArray, r0::Integer, c0::Integer, rn::Integer, cn::Integer, h0::Real )
        Δr = rn - r0
        Δc = cn - c0
        th0 = dtm[r0, c0] + h0
        steps = max( abs(Δr), abs(Δc) )
        r_inc = Δr / steps
        c_inc = Δc / steps
        heigths_profile = []
        coords_profile = []
        for (r, c) in zip( r0:r_inc:rn, c0:c_inc:cn ) 
            rint, cint = round.(Int64, [r, c])
            push!( heigths_profile, dtm[rint, cint][1] )
            push!( coords_profile, Tuple(GeoArrays.coords( dtm, [rint, cint])) )
            if dtm[rint, cint] >= th0
                break
            end
        end
        return heigths_profile, coords_profile
    end

    function generate_profiles( dtm::AbstractArray, x0::Real, y0::Real, h0::Real )
        values_profiles = []
        indexes_profiles = []
        r0, c0 = toIndexes(dtm, x0, y0)
        rm = size(dtm, 1)
        for α in [0, 90, 180, 270], β in 0:89
            rn, cn = rotate_point(rm, c0, r0, c0, α+β)
            values, indexes = getprofile(dtm, r0, c0, rn, cn, h0)
            push!(values_profiles, values)
            push!(indexes_profiles, indexes)
        end
        return values_profiles, indexes_profiles
    end

    function veiwshed( profile::AbstractVector )::AbstractVector
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
    end
=#





# === OLD VIEWSHED TEST ===========================================================================================================================

profile1 = [
    (0,12),
    (25,12),
    (50,8),
    (75,10),
    (100,13),
    (125,11),
    (150,15),
    (175,14),
    (200,23)
]
profile2 = [
    (0,500),
    (25,12),
    (50,8),
    (75,10),
    (100,9),
    (125,11),
    (150,15),
    (175,14),
    (200,23)
]
profile3 = [
    (0,12),
    (25,12),
    (50,8),
    (75,10),
    (100,16),
    (125,11),
    (150,500),
    (175,14),
    (200,23)
]
profile = profile3

visible = [ profile[1], profile[2] ]
vali = abs( ( profile[2][2] - profile[1][2] ) / ( profile[2][1] - profile[1][1] ) )
for i in 3:length(profile)
    println("vali: $vali")
    val = ( profile[i][2] - profile[i-1][2] ) / ( profile[i][1] - profile[i-1][1] )
    print("$(profile[i]): $val")
    if val >= abs(vali)
        print("    PUSH")
        push!( visible, profile[i] )
    end
    println("\n")
    vali = vali < 0 ? vali+val : vali-val
end
visible

# ============================================================================================================================================
