for i in 1:5
   print(i, ", ")
end
include("./src/coalescent.jl")

using PhyloNetworks
using PhyloPlots
using RCall
#net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
net = readTopology("((A:7,((B:2,C:2):3)#H1:2::0.4):3,(#H1:1::0.6,D:6):4);")
m=PhyloNetworks.mulTree(net)
preorder!(m)
a=PhyloNetworks.parentTrees(net,resetNodes=false)
a[3]

using Combinatorics
collect(partitions(a))

tops=unlabeledGenerate(4)
splitsHistories([4,4],tops)

b=descendenceMatrix(a[2])

w=[0.16,0.24,0.24,0.36];
w=[0.2,0.25,0.25,0.3]


w=[0.25,0.25,0.25,0.25]
w=[0.5,0.0,0.0,0.5]
df=PhyloNetworks.vcvParent(net,w)


mu=[0,0,0,0]
sig=convert(Matrix,df)

using Distributions
using DataFrames
mydist=MvNormal(mu,sig)
b=rand(mydist)

b=ParamsBM(3, 1)

using GLM
dat=DataFrame(b=b,tipNames=string.(names(df)))
fitTrait1 = PhyloNetworks.phyloParentlm(@formula(b ~ 1), dat,net,df)



PhyloNetworks.descendenceMatrix(a[2])


df=DataFrame(zeros(Float64,4,3))
names!(df,Symbol.(names))


PhyloNetworks.removeClade(net,net.node[5],true)

x=deepcopy(net)
mat=PhyloNetworks.unlabeledGenerate(2);
y=PhyloNetworks.decompILS!(x,x.node[4],1.0,mat)

PhyloNetworks.levelorder!(net)
plot(net, :R, showEdgeNumber=true, showNodeNumber=true)
plot(a[1], :R, showEdgeNumber=true, showNodeNumber=true)
plot(net, :R, showGamma=false,showNodeNumber=false,showEdgeLength=true);


"C:\\Users\\justison\\Desktop\\PhyloNetworks.jl"

PhyloNetworks.sharedPathMatrix(net)