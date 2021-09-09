module ValuableAdditions
#Module for functions thar are needed but not present in Julia

export rPaste #Functions

#=
Functions meant to imitate the behaviour of R's "paste"/"paste0" functions
=#

#Function for a single vector, returns aa vector of the elements converted to string 
function rPaste( a::AbstractVector, collapse::Union{String, Nothing} = nothing )
    if isnothing(collapse)
        return string.(a)
    else
        return join( a, collapse )
    end
end

#Function for two vectors, concatenates them element by element returning a vector of strings
 #if one of the vectors is greater that the other, the elements of the shorter one are repeated in
 #the concatenation
function rPaste( a::AbstractVector, b::AbstractVector, collapse::Union{String, Nothing} = nothing )
    if length(a) == length(b)
        res = @. string(a)*string(b)
        return isnothing(collapse) ? res : join( res, collapse )
    else
        la, lb = length(a), length(b)
        i, j = 1, 1
        if la > lb
            res = similar( a, String )
            while i <= la
                res[i] = string(a[i])*string(b[j])
                i += 1
                j = j+1 > lb ? 1 : j+1
            end
        else
            res = similar( b, String )
            while i <= lb
                res[i] = string(a[j])*string(b[i])
                i += 1
                j = j+1 > la ? 1 : j+1
            end
        end
        return isnothing(collapse) ? res : join( res, collapse )
    end
end

#Function that attaches String "b" to the end of each element of Vector "a", joining the result in a single
 #string with the String "sep" as separator    
function rPaste( a::AbstractVector, b::AbstractString, sep::AbstractString )
    arr = @. string( string(a), b )
    return join( arr, sep )
end



#Function to find all indexes of matching String or Regex
function findindexes( needle::Union{AbstractString, AbstractPattern, AbstractChar}, haystack::AbstractString )
    indexes = []
    val = findfirst( needle, haystack )
    while true        
        if isnothing(val) break end
        
        i = val[1]
        append!( indexes, i )
        i += 1
        val = findnext( needle, haystack, i )
    end
    return indexes
end

end # module