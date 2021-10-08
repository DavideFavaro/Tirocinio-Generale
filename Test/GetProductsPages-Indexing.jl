function count( maxNumber, count )
    if !isnothing( maxNumber ) && maxNumber > 1 && maxNumber < count
        count = maxNumber
    end
    return [ count÷100, count%100 ]
end


function indexing( count )
    if count[1] > 0
        for i in range( 0, count[1]-1, step=1 )
            println( "start=$(i*100) & rows=100" )
            println( "...\\$(i+1).xml")
        end
    end
    if count[2] > 0
        println( "start=$(count[1]*100) & rows=$(count[2])" )
        println( "...\\$(count[1]+1).xml" )
    end
end



function indexing2( start::Integer, val::Integer, maxNumber::Union{Nothing, Integer}=nothing )
    count = val
    println( "COUNT : ", count )

#Check if `start` has a sensible value
    if start > val
        start = val - 1
    end
    println( "START : ",start)

# Check if `maxNumber` has a sensible value
    if !isnothing( maxNumber ) && maxNumber > 0 && maxNumber < val
        count = maxNumber
    end

    if count > start
        count -= start
    end

    println( "COUNT(U1) : ", count )
    count = [ count ÷ 100, count % 100 ]
    println( "MAXNUMBER : ", maxNumber )
    println( "COUNT(U2) : ", count )
    println()

    if count[1] > 0
        println("count[1] > 0")
        for i in range( 0, count[1] - 1, step=1 )
            println( "start=$( start + ( i * 100 ) ) & rows=100" )
        end
    end
    if count[2] > 0
        println("count[2] > 0")
        println( "start=$(count[1] * 100 + start) & rows=$(count[2])" )
    end
    println()
    println("THE END")
end