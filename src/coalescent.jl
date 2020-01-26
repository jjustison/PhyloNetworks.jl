function g_ij(i::Int64,j::Int64,t)
    ##OverFlowError if i>20
    
    s=0 ##The sum of each term
    for k in j:i
        e=exp(-k*(k-1)*t/2)
        numerator=((2*k)-1)*((-1)^(k-j))
        denominator=factorial(j)*factorial(k-j)*(j+k-1)
    
        mult_term=1
        for y in 0:(k-1)
            mult_term*=(((j+y)*(i-y))/(i+y))
        end
        val=e*(numerator/denominator)*mult_term
        s+=val
    end
    return(s)
end

function nHistories(i::Int64,j::Int64)
    ##the number of histories where i lineages coalesce into j lineages 
    
    mult_term=1
    (j==1) && (j=2) ## coalescing into 1 or 2 lineages have the same number of histories but binomial(k,2) returns 0 when k=1
    for k in (j+1):i
        mult_term*=binomial(k,2)
    end
    return(mult_term)
end

function nTops(n::Int64,labeled=true)
    if labeled ##Number of labeled tree topologies
        n==1 && return(1)
        return(Int(factorial(2*n-3)/((2^(n-2)*factorial(n-2)))))
    
    else ##number of unlabeled topologies i.e. tree shapes
        n==1 && return(1) ##base case for recursion
        odd=true
        (mod(n,2)==0) && (odd=false)

        if odd ## n is odd
            val=0
            for p in 1:((n-1)/2)
                val+=nTops(p,false)*nTops(n-p,false)
            end
        else ## n is even
            t2=nTops(n/2,false)
            term=t2*(t2+1)/2

            s=0
            for p in 1:((n/2)-1)
                s+=nTops(p,false)*nTops(n-p,false)
            end
            val=term+s
        end
        return(Int(val))
    end
end

function nLabeledShape(n::Int64,s::Int64)
    ##Compute the number of labeled topologies of a given tree shapes
    ##Where:
        #n= # of taxa
        #s= # of symetric nodes i.e. nodes where its two subtrees have the same shape
    return(factorial(n)/(2^s))
end

function nHistoriesTop(desc_nodes::Array{Int64,1})
    return(factorial(length(desc_nodes))*prod((x->(1/(x+1)).(desc_nodes))))
end

function unlabeledGenerate(n)
    vect_tops=Dict{Int64, Array{Array{Int64,1},} }()
    sym_tops=Dict{Int64,Array{Int64,1}}()
    prod_tops=Dict{Int64,Array{Float64,1}}()

    ##Initialize the first two sets
    vect_tops[1]=[[]]
    sym_tops[1]=[0]
    prod_tops[1]=[1]

    vect_tops[2]=[[0]]
    sym_tops[2]=[1]
    prod_tops[2]=[1]

    for i in 3:n
        vect_tops[i]=[]
        sym_tops[i]=[]
        prod_tops[i]=[]

        for j in 1:(floor(i/2))
            ##make new shapes by combining all the the j_th and (n-j)_th arrays

            bottom=j
            top=i-j

            for k in 1:(length(vect_tops[bottom]))
                for l in 1:(length(vect_tops[top]))
                   first=vect_tops[bottom][k] 
                   second=vect_tops[top][l]

                   c=Int64[]
                   append!(c,first)
                   append!(c,second)
                   push!(c,length(c))
                   push!(vect_tops[i],c)

                   sym_new=(sym_tops[bottom][k]+sym_tops[top][l]) ##The number of symmetry points is the sum of the symmetry points of the two subtrees
                   (first==second) && (sym_new+=1) ## If the two shapes are the same then this makes another point of symmetry
                   push!(sym_tops[i],sym_new)

                   push!( prod_tops[i],prod( (x->(1/(x+1))).(c) ) ) ##TODO could make this more efficient by instead computing a new w from the c vector, just multiply the to previous products together
                end
            end
        end
    end
    return ((vect_tops,sym_tops,prod_tops))
end


function splitsHistories(splits::Array{Int64,1}, tops_mat)
    total=[0]

    spltsHisRec(splits,1,1,tops_mat,total)
   
    return total[1]

end

function spltsHisRec(splits::Array{Int64,1},ind_splits,running_val,tops_mat,total::Array{Int64,1} )
    split=splits[ind_splits]
    println("We're in split number ",ind_splits)
    println("the running value is: ", running_val)
    ##Get all the symmetries and products for the split we're interested in  
    syms=tops_mat[2][splits[ind_splits]]
    prods=tops_mat[3][splits[ind_splits]]

    nshapes=length(syms) ##the number of treeshapes for a given split

    for shape in 1:nshapes ##Do calculation for each shape
        val=running_val ##reset the value for each shape

        val*=factorial(split)/(2^(syms[shape]))
        val*=prods[shape]

        if ind_splits==length(splits) ##Base Case: we went thru each split and did the calculation
            val*=factorial(sum((x->x-1).(splits)))
            total[1]+=val 
            println("we're adding ", val, " to the total")
        
        else ##we still have more splits to do calculations
            spltsHisRec(splits,ind_splits+1,val,tops_mat,total)
        end
    end
end




