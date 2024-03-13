#pushfirst!(LOAD_PATH,"C:/Users/justison/Desktop/repos/PhyloNetworks.jl")
using Revise
using PhyloNetworks
using PhyloPlots
#using RCall
using Distributions
using DataFrames
using stats

using GLM
#using Plots
#using Combinatorics

net = readTopology("(O:4,(A:3,((B:0.4)#H1:1.6::0.92,((#H1:0::0.08,C:0.4):0.6,(D:.2,E:.2):0.8):1):1):1);");
net.edge[5].gamma

net= readTopology("(((A:3,((B:1,C:1):1)#H1:1::0.5):1,(#H1:1::0.5,D:3):1):1,((E:3,((F:1,G:1):1)#H2:1::0.5):1,(#H2:1::0.5,H:3):1):1);");
net= readTopology("((A:3,((B:1,C:1):1)#H1:1::0.5):1,(#H1:1::0.5,D:3):1);");
#preorder!(net)

pts=PhyloNetworks.parentTrees(net,report_hybsorting=true)
hyb_sortings = (x->x[2]).(pts)
weights=collect(zip(repeat([1/4],4),hyb_sortings)) #for larger network
weights=collect(zip([0.5,0,0,0.5],hyb_sortings))
#weights=collect(zip([0.36,0.24,0.24,0.16],hyb_sortings)) #equal to when gammas are 0.6,0.4
#weights=collect(zip([0.36,0.24,0.24,0.16],hyb_sortings)) # equal to when gammas are 0.5,0.5
#weights=collect(zip([0.36,0.24,0.24,0.16],hyb_sortings)) #equal to when gammas are 0.4,0.4
x=PhyloNetworks.vcvParent3(net,weights)
vcv(net)

plot(net,showedgenumber=true,shownodenumber=true)


v_mat = PhyloNetworks.sharedPathMatrix(net)


for nd in net.hybrid
    e=PhyloNetworks.getMajorParentEdge(nd)
    PhyloNetworks.setGamma!(e,0.5)
end

x = PhyloNetworks.mulTree(net)

scalar=0.0001
ngenes = 1:5:100
nreps = 10000

##get the theoretical pt probs
pt_probs = PhyloNetworks.parentTreeProbs(net;scalar=1/10)
w = (x->x[2]).(pt_probs)


sum((x-> x[2]).(pt_probs))

per_gene_ests= Array{Array{Float64,1},1}()
for ngene in ngenes
    ests=zeros(nreps);
    println(ngene)
    for i in 1:nreps
        println([ngene,i])
            ##simualte the trait 
        b,VCV,w_obs=traitPTSim(net,ngene,w,1);

        ##Fit Model
        dat=DataFrame(b=b,tipNames=string.(names(VCV)));
        fitTrait1 = phylolm(@formula(b ~ 1), dat, net); ##Fitting the traits on the Network - Wrong Model
        #fitTrait1 = PhyloNetworks.phyloParentlm(@formula(b ~ 1), dat, net,VCV); ##Fitting the traits on the parent trees - Correct Model
        sigma_est=sigma2_estim(fitTrait1); ##get estimates of the evolutionary rate
        ests[i]=sigma_est;
    end
    push!(per_gene_ests,ests);
end


vars = (x->var(x)).(per_gene_ests) 
means = (x->mean(x)).(per_gene_ests)
dist = (x->mean(abs.(x.-1))).(per_gene_ests)
Plots.plot(1:5:100,vars)

total=Float64[]
for i in per_gene_ests
    append!(total,i)
end
Plots.scatter(sort(repeat(1:100,1000)),total)


###############################
#### Test population Sizes ####
###############################


pop_sizes=tan.(range(atan(1/100),atan(100),length=30))
nreps = 10


per_gene_ests1= Array{Array{Float64,1},1}()
per_gene_ests2= Array{Array{Float64,1},1}()
for pop in pop_sizes
    wrong_ests=zeros(nreps);
    correct_ests=zeros(nreps);
    println(pop)
    for i in 1:nreps
        println([pop,i])

        ##get the theoretical pt probs
        pt_probs = PhyloNetworks.parentTreeProbs(net;scalar=pop)
        
        ##simualte the trait 
        b,VCV,w_obs=traitPTSim(net,Inf,pt_probs,1);
        V=Matrix(VCV)
        ##Fit Model
        dat=DataFrame(b=b,tipNames=string.(names(VCV)));
        V_orig = vcv(net)
        fitTrait1 = phylolm(@formula(b ~ 1), dat, net,reml=true) ##Fitting the traits on the Network - Wrong Model
        fitTrait2 = phylolm(@formula(b ~ 1), dat, net,model="BMVCV",VCV=Matrix(vcv(net))) ##Fitting the traits on the parent trees - Correct Model
        #sigma_est1=sigma2_estim(fitTrait1); ##get estimates of the evolutionary rate
        sigma_est2=sigma2_estim(fitTrait2); ##get estimates of the evolutionary rate
        #wrong_ests[i]=sigma_est1;
        #correct_ests[i]=sigma_est2;
    end
    #push!(per_gene_ests1,wrong_ests);
    #push!(per_gene_ests2,correct_ests);
end



vars = (x->var(x)).(per_gene_ests1) 
means = (x->mean(x)).(per_gene_ests1)
dist = (x->mean(abs.(x.-1))).(per_gene_ests1)
Plots.plot(1:30,means)








##simulate the trait
function traitPTSim(net,ngene,weights::Array{Tuple{Float64,Dict{Int64, Dict{Int64,Set{String}}}},1},sigma2=1)
    probs= (x->x[1]).(weights)
    hyb_sortings = (x->x[2]).(weights)

    if ngene==Inf
        w_obs=probs
    else
        ##set up the distribution for genes along parent trees
        gene_samp_dist=Multinomial(ngene,w) 
        ##sample genes for 'observed' proportions
        w_obs=rand(gene_samp_dist)/ngene
    end


    new_weights=collect(zip(w_obs,hyb_sortings))
    ##set up mvNorm for simulating traits
    mu=zeros(length(net.leaf))
    VCV=PhyloNetworks.vcvParent3(net,new_weights)
    sigma=sigma2.*Matrix(VCV)
    mvdist=MvNormal(mu,sigma)
    ##Simulate traits
    traits=rand(mvdist)

    return (traits,VCV,w_obs)
end


sigma2=3
mu=zeros(length(net.leaf))
V= vcv(net)
sigma=sigma2.*Matrix(V)
mvdist=MvNormal(mu,sigma)

CSV.write("test_V.csv",V)


nreps=10000
ests=[]
for rep in 1:nreps

    traits=rand(mvdist) #simulate traits
    dat=DataFrame(b=traits,tipNames=string.(names(V)));
    CSV.write("testdat.csv",dat)
    fitTrait2 = phylolm(@formula(b ~ 1), dat, net,model="BMVCV",VCV=Matrix(vcv(net))) #fit correct model
    sigma2_phylo(fitTrait2)
    push!(ests,sigma2_phylo(fitTrait2))
end