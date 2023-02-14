using PhyloNetworks
using Distributions
using DataFrames
using GLM

##the tree
net = readTopology("((A:7,((B:2,C:2):3)#H1:2::0.5):3,(#H1:1::0.5,D:6):4);")
ngene=30
##high Ne case 
#w=[0.25,0.25,0.25,0.25]

##low Ne case
w=[0.5,0.0,0.0,0.5]

reps=100000
ests=zeros(reps)
for i in 1:reps
    
    ##simualte the trait 
    b,VCV=traitPTSim(net,ngene,w,1)

    ##Fit Model
    dat=DataFrame(b=b,tipNames=string.(names(VCV)))
    fitTrait1 = phyloNetworklm(@formula(b ~ 1), dat, net) ##Fitting the traits on the Network - Wrong Model
    #fitTrait1 = Phylonetworks.phyloParentlm(@formula(b ~ 1), dat, net,VCV) ##Fitting the traits on the parent trees - Correct Model
    sigma_est=sigma2_estim(fitTrait1) ##get estimates of the evolutionary rate
    ests[i]=1-sigma_est
end












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