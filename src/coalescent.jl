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



function nTops(n::Int64,labeled=true::Bool,rooted=true::Bool)
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
    ##we will create dictionaries to store the bits we want for unlabeled tree shapes
    ##Each key in the dictionary will point to a vector of the bits for each tree shape 
        ##e.g. sym_tops[5] will give us a vector of how many points of symmetry each tree shape with 5 taxa has 
    vect_tops=Dict{Int64, Array{Array{Int64,1},} }() ##this bit is the node vector representation of the tree shape. Each element in the vector represents an internal node and it denotes how many internal nodes it has as descendents
    sym_tops=Dict{Int64,Array{Int64,1}}() ##this bit stores the number of reflection points for a given tree shape
    prod_tops=Dict{Int64,Array{Float64,1}}() ##if x is the node vector representation of the tree shape, this bit stores the product of 1/(1+x)

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
    return ((vect_tops,sym_tops,prod_tops));
end

function splitsHistories(splits::Array{Int64,1}, tops_mat)
    total=[0]

    spltsHisRec(splits,1,1,tops_mat,total)
   
    return total[1];

end

function spltsHisRec(splits::Array{Int64,1},ind_splits,running_val,tops_mat,total::Array{Int64,1} )
    split=splits[ind_splits]
    ##Get all the symmetries and products for the split we're interested in  
    syms=tops_mat[2][splits[ind_splits]]
    prods=tops_mat[3][splits[ind_splits]]

    nshapes=length(syms) ##the number of treeshapes for a given number of taxa in the split

    for shape in 1:nshapes ##Do calculation for each shape
        val=running_val ##reset the value for each shape

        val*=factorial(split)/(2^(syms[shape]))
        val*=prods[shape]

        if ind_splits==length(splits) ##Base Case: we went thru each split and did the calculation
            val*=factorial(sum((x->x-1).(splits)))
            total[1]+=val 
        else ##we still have more splits to do calculations
            spltsHisRec(splits,ind_splits+1,val,tops_mat,total)
        end
    end
end

function findLowNode(net::HybridNetwork) ##TODO. Make this generate a list of all the nodes to do in order instead of one at a time. Also, since net.node never changes until the end, we could have it store the node indices instead of pointers to the nodes themselves. 
    numbers=(x->x.number).(net.node) ##store old numbers '
    nd_changed= net.nodes_changed ##store node ordering
    rankNodes!(net)
    

    
    ##get all hybrid nodes
    hybs=net.hybrid

    nds=Set() ##This will have all the nodes we need to decompose
    for hyb in hybs
        union!(nds,Set(preorder(net,hyb))) ##add the hyb node and all its descendents to the list
    end
    nds=collect(nds) ##turn the nodes into an array
    nds= nds[(x->(!x.leaf)).(nds)] ##get rid of tips
     
    nds=nds[sortperm((x->x.number).(nds))]
    ((x,y)->x.number=y).(net.node,numbers) ##reset node numbers to their old values
    net.nodes_changed=nd_changed ##reset node ordering


    return nds
end

function decompILS(network::HybridNetwork,decomp_node::Node,prb::Float64,tm)

    ## Make a copy of the network that we can muck around with.
    nd_ind=findfirst(x->x==decomp_node,network.node)
    net=deepcopy(network)
    nd=net.node[nd_ind]

    par=getParents(nd)[1] 
    nmes=(x->x.name).(getChildren(nd))

    c_edge=getChildrenEdges(nd)[1].length   ##get the edge length of the first child, ignore other child edge lengths, assuming ultrametricy they should all be the same anyways
    p_edge=getChildrenEdges(par)[1].length  ##get the edge length of the parent
    elen=c_edge+p_edge ##This will be the edge length we assign to our coalesced edges

    removeClade!(net,nd,false)
    par_ind=findfirst(x->x==par,net.node)

    subsets=collect(partitions(nmes)) ##we need the power set (excluding the empty set) of all names for decomposition 

    ##We ultimately want ot report the decomposed trees and their respective probabilities after everything is all said and done
    decomposed_trees=HybridNetwork[]
    probs=Float64[]




    for subset in subsets
        ##start with a fresh tree for adding things back
        decomp_net=deepcopy(net)
        decomp_par=decomp_net.node[par_ind]

        new_node_num=maximum((x-> x.number).(net.node))
        new_edge_num=maximum((x-> x.number).(net.edge))

        nms_coal=String[] ##These are the names of the coalesced tips we're adding to the tree 
        for nameset  in subset ##We need the individual sets within subset to concatenate the names 
            if length(nameset)==1
                nm=nameset[1]
            else
                nm=join(nameset,"|")
            end
            push!(nms_coal,nm)
        end

        for nm in nms_coal #Sequentially add each coalesced tip back to the tree
            ##create a new node and edge for the tip
            ntip=Node(new_node_num,true)
            ntip.name=nm
            nedge=Edge(new_edge_num,elen)
            new_edge_num+=1
            new_node_num+=1

           
            ##make appropriate connections
            push!(decomp_par.edge,nedge)
            push!(nedge.node,ntip)
            push!(nedge.node,decomp_par)
            push!(ntip.edge,nedge)

            ##add the new nodes and edges to the network
            push!(decomp_net.node,ntip)
            push!(decomp_net.edge,nedge)
            push!(decomp_net.leaf,ntip)
        end 
        
        ##update the attributes of the network
        num_added=length(nms_coal)
        decomp_net.numTaxa+=num_added
        decomp_net.numNodes+=num_added
        decomp_net.numEdges+=num_added
        preorder!(decomp_net) ##preorder things for good measure 

        ##Compute the probability of the tree
        gij=g_ij(length(nmes),length(nms_coal),p_edge) ## the probability of i lineages coalescing into j lineages by time t
        hp=splitsHistories(length.(subset),tm) ##The number of coalescent histories that lead to these partitions
        ch=nHistories(length(nmes),length(nms_coal)) ##The total number of coalescent histories for i coalescing into j lineages

        treeprob=gij*(hp/ch)*prb ##TODO I'm concerned with underflow so make this deal with small numbers. perhaps switch to log probs
        treeprob==0 && error("Likely underflow of the probability")
        
        push!(decomposed_trees,decomp_net)
        push!(probs,treeprob)
    end
    return collect(zip(decomposed_trees,probs))
end

function emptyHybSorting(net::HybridNetwork;type="name")
    if type == "name"
        hybsort= Dict{Int64, Dict{Int64,Set{String}}}()
        hybs=net.hybrid
        for hyb in hybs
            hybsort[hyb.number]=Dict{Int64,Set{String}}()
            for nd in getParents(hyb)
            hybsort[hyb.number][nd.number]=Set{String}() 
            end
        end
    elseif type == "node"
        hybsort= Dict{Int64, Dict{Int64,Set{Node}}}()
        hybs=net.hybrid
        for hyb in hybs
            hybsort[hyb.number]=Dict{Int64,Set{Node}}()
            for nd in getParents(hyb)
            hybsort[hyb.number][nd.number]=Set{Node}() 
            end
        end
    else
        error("we didn't get a valid type. Either 'node' or 'name' is accepted")
    end
    return hybsort
end

function decompHyb(network::HybridNetwork,nd::Node,prb::Float64,hyb_sorting::Dict{Int64, Dict{Int64,Set{String}}} )
    ## Make a copy of the network that we can muck around with.
    nd_ind=findfirst(x->x==nd,network.node)
    net=deepcopy(network)
    nd=net.node[nd_ind]
    
    hyb_num= nd.number

    lpar=getMajorParent(nd) ##left parent
    rpar=getMinorParent(nd) ##right parent
    lnum=lpar.number 
    rnum=rpar.number
    lp_e=getMajorParentEdge(nd) 
    rp_e=getMinorParentEdge(nd)  

    c_elen=getChildrenEdges(nd)[1].length   ##get the edge length of the first child, ignore other child edge lengths, they should be the same anyway
    lp_elen=lp_e.length
    rp_elen=rp_e.length

    l_gamma=lp_e.gamma
    r_gamma=rp_e.gamma

    removeClade!(net,nd,false)

    lpar_ind=findfirst(x->x==lpar,net.node) ##The index of the left parent in net.node
    rpar_ind=findfirst(x->x==rpar,net.node) ##The index of the right parent in net.node

    ##create all combinations for how the lineages could sort between the two hybrid getParents
    nmes=(x->x.name).(getChildren(nd)) ##get each lineage
    lsort=collect(combinations(nmes)) ##Generate all possible ways that things could sort to the left
    push!(lsort,[]) ## the function above doesn't account for the fact that no species could sort left
    rsort=(x->setdiff(nmes,x)).(lsort) ## everything that doesn't sort left must go right. NOTE: we assume that no lineages have the same name here
    combs=collect(zip(lsort,rsort)) ##associate the proper sortings together to generate all ways that lineages can sort

    ##These will store the decomposed trees and their probabilities 
    decomposed_trees=HybridNetwork[]
    probs=Float64[]
    hybsorts=Dict{Int64, Dict{Int64,Set{String}}}[]

    for i in combs
        ##start with a fresh slate for adding things
        hybsort=deepcopy(hyb_sorting)
        decomp_net=deepcopy(net)
        ldecomp_par=decomp_net.node[lpar_ind]
        rdecomp_par=decomp_net.node[rpar_ind]
        
        leftgoing=first(i)
        rightgoing=last(i)


        ####################
        ###Deal with left###
        ####################
        if length(leftgoing)>1 ##if more than 1 tip, we need to make a node for the mrca of all tips going left
            nedge=Edge(-1,lp_elen)
            ##Make appropriate connections
            push!(ldecomp_par.edge,nedge)
            push!(nedge.node,ldecomp_par)
            ldecomp_par=Node(-1,false,false)
            push!(ldecomp_par.edge,nedge)
            pushfirst!(nedge.node,ldecomp_par)
            ##add the new nodes and tips to the network
            push!(decomp_net.node,ldecomp_par)
            push!(decomp_net.edge,nedge)
            ##update attributes of the network
            decomp_net.numNodes+=1
            decomp_net.numEdges+=1
            r_elen=c_elen
        elseif  length(leftgoing)==0 ##if there is nothing going to the left we need to resolve the parent node 
            if ldecomp_par.hybrid
                resolveLonelyHybrid!(decomp_net,ldecomp_par,true)
            else
                fuseedgesat!(findfirst(x->x==ldecomp_par,decomp_net.node),decomp_net)
            end
        else ##only 1 tip going left
            r_elen=c_elen+rp_elen
        end
        for nm in leftgoing ##add tips going left to the network
            ntip=Node(-1,true)
            ntip.name=nm
            nedge=Edge(-1,lp_elen+c_elen)

             ##make appropriate connections
             push!(ldecomp_par.edge,nedge)
             push!(ntip.edge,nedge)
             push!(nedge.node,ntip)
             push!(nedge.node,ldecomp_par)

            ##add the new nodes and tips to the network
            push!(decomp_net.node,ntip)
            push!(decomp_net.edge,nedge)
            push!(decomp_net.leaf,ntip)
        end

        ####################
        ###Deal with right##
        ####################
        if length(rightgoing)>1 ##if more than 1 tip, we need to make a node for the mrca of all tips going right
            nedge=Edge(-1,rp_elen)
            ##Make appropriate connections
            push!(rdecomp_par.edge,nedge)
            push!(nedge.node,rdecomp_par)
            rdecomp_par=Node(-1,false,false)
            push!(rdecomp_par.edge,nedge)
            pushfirst!(nedge.node,rdecomp_par)
            ##add the new nodes and tips to the network
            push!(decomp_net.node,rdecomp_par)
            push!(decomp_net.edge,nedge)
            ##update attributes of the network
            decomp_net.numNodes+=1
            decomp_net.numEdges+=1
            r_elen=c_elen
        elseif  length(rightgoing)==0 ##if there is nothing going to the right we need to resolve the parent node 
            if rdecomp_par.hybrid
                resolveLonelyHybrid!(decomp_net,rdecomp_par,true)
            else
                fuseedgesat!(findfirst(x->x==rdecomp_par,decomp_net.node),decomp_net)
            end
        else #only 1 tip going right
            r_elen=c_elen+rp_elen
        end
        for nm in rightgoing ##add tips to the right parent
            ntip=Node(-1,true)
            ntip.name=nm
            nedge=Edge(-1,r_elen)

             ##make appropriate connections
             push!(rdecomp_par.edge,nedge)
             push!(ntip.edge,nedge)
             push!(nedge.node,ntip)
             push!(nedge.node,rdecomp_par)

            ##add the new nodes and tips to the network
            push!(decomp_net.node,ntip)
            push!(decomp_net.edge,nedge)
            push!(decomp_net.leaf,ntip)
        end
        ##the tree should now be deocmposed 
        
        ##Update the attributes of the network
        num_added=length(rightgoing)+length(leftgoing)
        decomp_net.numNodes+=num_added
        decomp_net.numEdges+=num_added
        decomp_net.numTaxa+=num_added
        preorder!(decomp_net) ##preorder things for good measure


        ##record which lineages went left and right along the reticulation
        left_names=String[]
        right_names=String[]
        lnms= (x->split(x,"|")).(leftgoing)
        rnms= (x->split(x,"|")).(rightgoing)
        println(lnms)
        println(rnms)
        for nms in lnms
            append!(left_names,String.(nms))
        end
        for nms in rnms
            append!(right_names,String.(nms))
        end
        hybsort[hyb_num][lnum]=Set(left_names)
        hybsort[hyb_num][rnum]=Set(right_names)

        ##Compute the probability of the tree
        treeprob=(r_gamma^length(rightgoing))*(l_gamma^length(leftgoing))*prb ##TODO I'm concerned with underflow so make this deal with small numbers. perhaps switch to log probs
        treeprob==0 && error("Likely underflow of the probability")

        push!(decomposed_trees,decomp_net)
        push!(probs,treeprob)
        push!(hybsorts,hybsort)
    end
    return collect(zip(decomposed_trees,probs,hybsorts))
end

