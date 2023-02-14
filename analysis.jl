function matDist(A,B) ##compute the distance between two matrices
    ##Distances are compuated as the Frobenius norm of A-B, ||A-B||  
    norm(A-B,2)
end

##A lineage pt is pt1 and 0.4 
##D lineage pt is pt4 and 0.6
function parentTreeProbs4(t,gmma)
    pt1=(1-exp(-t))*(gmma)+exp(-t)*(gmma^2)
    pt2=exp(-t)*gmma*(1-gmma)
    pt3=exp(-t)*gmma*(1-gmma)
    pt4=(1-exp(-t))*((1-gmma))+(exp(-t)*((1-gmma)^2))

    return [pt1,pt2,pt3,pt4]
end

function ptp4simple(t,gmma)
    pt1=(1-t)*(gmma)+t*(gmma^2)
    pt2=t*gmma*(1-gmma)
    pt3=t*gmma*(1-gmma)
    pt4=(1-t)*((1-gmma))+(t*((1-gmma)^2))

    return [pt1,pt2,pt3,pt4]
end





V=Matrix(vcv(net))
dists=Float64[]

for val in range(0,100000,length=101)
    probs=parentTreeProbs4(val,0.4)

    Matrix(PhyloNetworks.vcvParent(net,probs))
    d=matDist(V,V2)
    push!(dists,d)
end


