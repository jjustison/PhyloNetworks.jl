pushfirst!(LOAD_PATH,"C:/Users/justison/Desktop/repos/PhyloNetworks.jl")
using Revise
using PhyloNetworks
using PhyloPlots
using RCall
using Distributions
using DataFrames
using GLM
using Plots

#net= readTopology("(((#H16:0.05368826116,#H27:0.3207744925):0.1006108928,((((t7:0.1305863974)#H23:0.03502822662,(t14:0.1449360234,(t9:0.1181312427)#H27:0.02680478076):0.02067860053):0.2459876316,((((t5:0.07857989615,(t10:0.05043758827,t8:0.05043758827):0.02814230788):0.1479715226,#H20:0.0503032565):0.03824463648)#H17:0.1204214189)#H16:0.02638478156):0.0745330051,((t1:0.1486468951,t13:0.1486468951):0.0110575952,#H23:0.02911809292):0.3264307704):0.05338136728):0.460483372,((((t6:0.1305965702)#H21:0.04565159199)#H20:0.1346353075,#H17:0.04608741449):0.1424901256,#H21:0.3227770251):0.5466264047);")
net = readTopology("((A:7,((B:2,C:2):3)#H1:2::0.5):3,(#H1:1::0.5,D:6):4);")
net = readTopology("C:/Users/justison/Documents/mynets/net29.txt")
net = readTopology("(((t19:0.2871026308)#H26:0.4947943372,(((((t1:0.09654463128,t4:0.09654463128):0.2531353214,#H28:0.0497256181):0.05737946088,#H26:0.05453036529):0.1114482259,((((t5:0.0004167988422,t15:0.0004167988422):0.199117527,t18:0.1995343258):0.1330880931,(((t10:0.1329623641,t20:0.1329623641):0.04526931012,t11:0.00364736851):0.1068498567,t8:0.0909196431):0.04754088798):0.07309595057,#H29:0.03324616104):0.11278927):0.1582198116,(((((t2:0.02147251407,t22:0.06553262051):0.1504194273,t12:0.2159520478):0.06385698246)#H39:0.0654570322,((t9:0.1329599145,(t23:0.06829953659)#H49:0.005606300993):0.08460752085,t24:0.2175674354):0.1276986271):0.2454392287,t3:0.5907052912):0.08602215993):0.1705959344):0.1526766145,((t25:0.1492511066,t17:0.2947087114):0.03053120729,((#H49:0.172600721)#H28:0.4202917547,((t14:0.0004755690989,#H39:0.04430258142):0.04836059676)#H29:0.3477738809):0.2133367394):0.0664171713);")
net = readTopology("(((t12:0.0406134002,((((t24:0.2744055017,t7:0.2744055017):0.2261360787,(t1:0.1015884872)#H42:0.1052548302):0.006580124502,t4:0.5071217049):0.01236765639,((t19:0.3213959723,#H44:0.1332227491):0.142121733)#H39:0.055971656):0.278516031):0.01423330077,(t21:0.4640005553,(t5:0.6385706426,t14:0.1115833917):0.07978742786):0.0938806226):0.187761307,(((t8:0.09887219626,t6:0.4411707113):0.1751632675,((t11:0.04741979486)#H45:0.3599017768,#H42:0.01203482139):0.2090124071):0.2716536659,((((((t20:0.05753635394,t10:0.05753635394):0.3278354229,#H45:0.337951982):0.2650128997,((t13:0.011916546)#H47:0.03313950524,(t22:0.03338578103,(t2:0.01919052856,#H47:0.007273982559):0.01419525247):0.01167027022):0.6053286252):0.09861838154,((t16:0.1881732232)#H44:0.2492610473,t9:0.4374342705):0.3115687875):0.08254943368,#H39:0.3680347864):0.01338700973,t3:0.03894281828):0.04304814313):0.1120123554);")
preorder!(net)

b = PhyloNetworks.blobInfo(net)

##Set the gamma values because I'm a lazy POS and haven't implemented writing inheritance probabilities in R yet
for nd in net.hybrid
    e=PhyloNetworks.getMajorParentEdge(nd)
    PhyloNetworks.setGamma!(e,0.5)
end



scalar=0.0001
ngenes = 1:100
nreps = 200

##get the theoretical pt probs
pt_probs = PhyloNetworks.parentTreeProbs(net;scalar=scalar)
w = (x->x[2]).(pt_probs)


sum((x-> x[2]).(pt_probs))

per_gene_ests= Array{Array{Float64,1},1}()
for ngene in ngenes
    ests=zeros(nreps)
    println(ngene)
    for i in 1:nreps
            ##simualte the trait 
        b,VCV,w_obs=traitPTSim(net,ngene,w,1)

        ##Fit Model
        dat=DataFrame(b=b,tipNames=string.(names(VCV)))
        fitTrait1 = phyloNetworklm(@formula(b ~ 1), dat, net) ##Fitting the traits on the Network - Wrong Model
        #fitTrait1 = Phylonetworks.phyloParentlm(@formula(b ~ 1), dat, net,VCV) ##Fitting the traits on the parent trees - Correct Model
        sigma_est=sigma2_estim(fitTrait1) ##get estimates of the evolutionary rate
        ests[i]=sigma_est
    end
    push!(per_gene_ests,ests)

end


vars = (x->var(x)).(per_gene_ests) 
means = (x->mean(x)).(per_gene_ests)
dist = (x->mean(abs.(x.-1))).(per_gene_ests)
Plots.plot(1:100,vars)

total=Float64[]
for i in per_gene_ests
    append!(total,i)
end
Plots.scatter(sort(repeat(1:100,1000)),total)

##simulate the trait
function traitPTSim(net,ngene,w,sigma2=1)

    if ngene==Inf
        w_obs=w
    else
        ##set up the distribution for genes along parent trees
        gene_samp_dist=Multinomial(ngene,w) 
        ##sample genes for 'observed' proportions
        w_obs=rand(gene_samp_dist)/ngene
    end

    ##set up mvNorm for simulating traits
    mu=zeros(length(net.leaf))
    VCV=PhyloNetworks.vcvParent(net,w_obs)
    sigma=sigma2.*convert(Matrix,VCV)
    mvdist=MvNormal(mu,sigma)
    ##Simulate traits
    traits=rand(mvdist)

    return (traits,VCV,w_obs)
end