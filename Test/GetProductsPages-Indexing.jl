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