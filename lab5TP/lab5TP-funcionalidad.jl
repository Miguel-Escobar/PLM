function genData( a, d )
    sz = zeros(Integer,nI+nJ)
    m = 1
    for i in I,
        sz[m] = size(d[i],1); m = m+1
    end
    for j in J
        sz[m] = size(a[j],1); m = m+1
    end
    
    nS = prod(sz)
    dmnd = Matrix{Integer}(undef,nI,nS)
    avlb = Matrix{Integer}(undef,nJ,nS)
    prob = ones(Float64,nS)
    
    s = 1
    itm = copy(sz)
    while s <= nS
        for i in I
            dmnd[i,s] = d[i][itm[i]]
            prob[s] = prob[s]*p[i][itm[i]]
        end
        for j in J
            avlb[j,s] = a[j][itm[j+nI]]
            prob[s] = prob[s]*q[j][itm[j+nI]]
        end    
    
        m = 1
        while m < nI+nJ
            if itm[m] > 1
                break
            end
            m = m+1
        end
        for w = 1:(m-1)
            itm[w] = sz[w]
        end
        itm[m] = itm[m]-1
        s = s+1
    end

    return dmnd, avlb, prob
end


function muestra( mdl )
    # Recibe un modelo que ha sido resuelto e imprime su status, valor objetivo, y valor de sus variables.
    println( mdl )
    st = termination_status( mdl )
    println( st )
    if st == MOI.OPTIMAL
        println( objective_value( mdl ) )
    end
    varArray = all_variables(mdl)
    for i = 1:size(varArray,1)
        println( varArray[i], " = ", value(varArray[i]) )
    end
    println()
end