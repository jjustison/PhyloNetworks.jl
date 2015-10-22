# functions written in classes.jl and moved here after tested
# for pseudolikelihood implementation (Stage2)
# Claudia August 2014
#
# in julia: include("functions.jl")

# tests of functions in examples_classes.jl outside git_laptop

# needed modules:
#using DataStructures # for updateInCycle with queue
using Base.Collections # for updateInCycle with priority queue
using DataFrames # for rep function and read/write csv tables
using NLopt # for branch lengths optimization
using GraphViz #for visualization

const DEBUG = false
const DEBUGC = false #more detailed prints
const CHECKNET = false #for debugging only
const REDIRECT = false # for debugging to a file later


include("auxiliary2.jl") #auxiliary2 for julia 0.4

include("update.jl")

include("undo.jl")

include("addHybrid.jl")

include("deleteHybrid.jl")

include("moves.jl")

include("readwrite.jl")

include("readData2.jl") #readData2.jl to try it in julia 0.4

include("optimization.jl")

include("pseudolik.jl")

include("visualization.jl")

include("descriptive.jl")

include("bootstrap.jl")
