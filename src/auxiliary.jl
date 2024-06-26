# auxiliary functions for all the other methods
# originally in functions.jl
# Claudia February 2015
#####################

function setCHECKNET(b::Bool)
    global CHECKNET
    CHECKNET = b
    CHECKNET && @warn "PhyloNetworks.CHECKNET is true: will slow snaq! down."
    b || @info "PhyloNetworks.CHECKNET set to false"
end

# ----- aux general functions ---------------

#based in coupon's collector: E+sqrt(V)
function coupon(n::Number)
    return n*log(n) + n
end

function binom(n::Number,k::Number)
    n >= k || return 0
    n == 1 && return 1
    k == 0 && return 1
    binom(n-1,k-1) + binom(n-1,k) #recursive call
end

function approxEq(a::Number,b::Number,absTol::Number,relTol::Number)
    if(a<eps() || b<eps())
        abs(a-b) < absTol
    else
        abs(a-b) < relTol*eps(abs(a)+abs(b))
    end
end

approxEq(a::Number,b::Number) = approxEq(a,b,1e-5,100)

# isEqual functions: to test if 2 edges (or 2 nodes etc.) "look" alike.
#                    Useful after a deepcopy of a network.
# For nodes (or edges etc.) in the same network, use instead n1 == n2 or n1 != n2.
function isEqual(n1::Node,n2::Node)
    return (n1.number == n2.number && approxEq(n1.gammaz,n2.gammaz) && n1.inCycle == n2.inCycle)
end

function isEqual(n1::Edge,n2::Edge)
    return (n1.number == n2.number && approxEq(n1.length,n2.length))
end

function isEqual(net1::HybridNetwork, net2::HybridNetwork)
    result = true
    result &= (net1.numTaxa == net2.numTaxa)
    result &= (net1.numNodes == net2.numNodes)
    result &= (net1.numEdges == net2.numEdges)
    ## result &= (net1.node == net2.node)
    ## result &= (net1.edge == net2.edge)
    result &= (net1.root == net2.root)
    result &= (net1.names == net2.names)
##    result &= (net1.hybrid == net2.hybrid)
    result &= (net1.numHybrids == net2.numHybrids)
##    result &= (net1.leaf == net2.leaf)
    result &= (net1.ht == net2.ht)
    result &= (net1.numht == net2.numht)
    result &= (net1.numBad == net2.numBad)
    result &= (net1.hasVeryBadTriangle == net2.hasVeryBadTriangle)
    result &= (net1.index == net2.index)
    result &= (net1.loglik == net2.loglik)
    return result
end


#------------- functions to allow for ------------#
#              missing values (lengths or gammas) #

# adds x+y but interprets -1.0 as missing: so -1.0 + x = -1.0 here.
function addBL(x::Number,y::Number)
    (x==-1.0 || y==-1.0) ? -1.0 : x+y
end
function multiplygammas(x::Number,y::Number)
    (x==-1.0 || y==-1.0) ? -1.0 : x * y
end

#------------- EDGE functions --------------------#

# warning: node needs to be defined as hybrid before adding to a
#          hybrid edge. First, an edge is defined as hybrid, and then
#          the nodes are added to it. If the node added is leaf, the
#          edge length is set unidentifiable (as it is external edge)
function setNode!(edge::Edge, node::Node)
    size(edge.node,1)  !=  2 || error("vector of nodes already has 2 values");
    push!(edge.node,node);
    if size(edge.node,1) == 1
        if edge.hybrid
            edge.isChild1 = node.hybrid
        end
        edge.istIdentifiable = !node.leaf
    else
        if node.leaf
            !edge.node[1].leaf || error("edge $(edge.number) has two leaves")
            edge.istIdentifiable = false;
        else
          if edge.hybrid
            if node.hybrid
                # @debug (edge.node[1].hybrid ? "hybrid edge $(edge.number) has two hybrid nodes" : "")
                edge.isChild1 = false;
	        else
	            edge.node[1].hybrid || error("hybrid edge $(edge.number) has no hybrid nodes");
	            edge.isChild1 = true;
	        end
          else #edge is tree
            if !edge.node[1].leaf
                if !node.hybrid && !edge.node[1].hybrid
                    edge.istIdentifiable = !edge.fromBadDiamondI
                else
                    if node.hybrid && (node.isBadDiamondI || node.isBadDiamondII || node.isBadTriangle)
                        edge.istIdentifiable = false
                    elseif edge.node[1].hybrid && (edge.node[1].isBadDiamondI ||edge.node[1].isBadDiamondII || edge.node[1].isBadTriangle)
                        edge.istIdentifiable = false
                    else
                        edge.istIdentifiable = true
                    end
                end
            else
                edge.istIdentifiable = false
            end
          end
        end
    end
end

# warning: node needs to be defined as hybrid before adding to a hybrid edge.
#          First, an edge is defined as hybrid, and then the nodes are added to it.
#          If there is a leaf in node, the edge.istIdentifiable=false
function setNode!(edge::Edge,node::Array{Node,1})
    size(node,1) ==  2 || error("vector of nodes must have exactly 2 values")
    edge.node = node;
    if(edge.hybrid)
      if(node[1].hybrid)
          edge.isChild1 = true;
      else
          node[2].hybrid || error("hybrid edge without hybrid node");
          edge.isChild1 = false;
      end
    end
    if(edge.node[1].leaf || edge.node[2].leaf)
        edge.istIdentifiable = false;
    else
        edge.istIdentifiable = true;
    end
end

"""
    getroot(net)

Node used to root `net`. If `net` is to be considered as semi-directed or
unrooted, this root node is used to write the networks' Newick parenthetical
description or for network traversals.

See also: [`isrootof`](@ref)
"""
getroot(net::HybridNetwork) = net.node[net.root]

"""
    isrootof(node, net)

`true` if `node` is the root of `net` (or used as such for network traversals
in case the network is considered as semi-directed); `false` otherwise.

    isleaf(node)
    isexternal(edge)

`true` if `node` is a leaf or `edge` is adjacent to a leaf, `false` otherwise.

See also: [`getroot`](@ref),
[`getparent`](@ref), [`getchild`](@ref)
"""
isrootof(node::Node, net::HybridNetwork) = node === getroot(net)

@doc (@doc isrootof) isleaf
isleaf(node::Node) = node.leaf

@doc (@doc isrootof) isexternal
isexternal(edge::Edge) = any(isleaf.(edge.node))

"""
    isparentof(node, edge)
    ischildof(node, edge)

`true` if `node` is the tail / head, or parent / child, of `edge`; `false` otherwise.
Assumes that the edge's direction is correct, meaning it's field `isChild1` is
reliable (in sync with the rooting).

See also: [`getparent`](@ref), [`getchild`](@ref), [`isrootof`](@ref)
"""
isparentof(node::Node, edge::Edge) = node === getparent(edge)
@doc (@doc isparentof) ischildof
ischildof( node::Node, edge::Edge) = node === getchild(edge)

"""
    hassinglechild(node)

`true` if `node` has a single child, based on the edges' `isChild1` field;
`false` otherwise.

See also: [`getchild`](@ref), [`getparent`](@ref)
"""
hassinglechild(node::Node) = sum(e -> getparent(e) === node, node.edge) == 1

"""
    getchild(edge)
    getchild(node)
    getchildren(node)

Get child(ren) **node(s)**.
- `getchild`: single child node of `edge`, or of `node` after checking that
  `node` has a single child.
- `getchildren`: vector of all children *nodes* of `node`.


    getchildedge(node)

Single child **edge** of `node`. Checks that it's a single child.

*Warning*: these functions rely on correct edge direction, via their `isChild1` field.

See also:
[`getparent`](@ref),
[`getpartneredge`](@ref),
[`isparentof`](@ref),
[`hassinglechild`](@ref).
"""
getchild(edge::Edge) = edge.node[edge.isChild1 ? 1 : 2]
getchild(node::Node) = getchild(getchildedge(node))

@doc (@doc getchild) getchildren
function getchildren(node::Node)
    children = Node[]
    for e in node.edge
        if isparentof(node, e)
            push!(children, getchild(e))
        end
    end
    return children
end

@doc (@doc getchild) getchildedge
function getchildedge(node::Node)
    ce_ind = findall(e -> isparentof(node, e), node.edge)
    length(ce_ind) == 1 || error("node number $(node.number) has $(length(ce_ind)) children instead of 1 child")
    return node.edge[ce_ind[1]]
end

function getchildrenedges(node::Node)
    child_edges=Edge[]
    for e in node.edge
        if node==getparent(e)
            push!(child_edges,e)
        end
    end
    return child_edges
end



"""
    getparent(edge)
    getparent(node)
    getparentminor(node)
    getparents(node)

Get parental **node(s)**.
- `getparent`: **major** (or only) parent node of `edge` or `node`
- `getparentminor`: minor parent node of `node`
- `getparents`: vector of all parent nodes of `node`.


    getparentedge(node)
    getparentedgeminor(node)

Get one parental **edge** of a `node`.
- `getparentedge`: major parent edge. For a tree node, it's its only parent edge.
- `getparentedgeminor`: minor parent edge, if `node` is hybrid
  (with an error if `node` has no minor parent).
If `node` has multiple major (resp. minor) parent edges, the first one would be
returned without any warning or error.

*Warning*: these functions use the field `isChild1` of edges.

See also: [`getchild`](@ref),
[`getpartneredge`](@ref).
"""
getparent(edge::Edge) = edge.node[edge.isChild1 ? 2 : 1]
@inline function getparent(node::Node)
    for e in node.edge
        if e.isMajor && ischildof(node, e)
            return getparent(e)
        end
    end
    error("could not find major parent of node $(node.number)")
end

@doc (@doc getparent) getparentminor
@inline function getparentminor(node::Node)
    for e in node.edge
        if !e.isMajor && node == getchild(e)
            return getparent(e)
        end
    end
    error("could not find minor parent of node $(node.number)")
end

@doc (@doc getparent) getparents
@inline function getparents(node::Node)
    parents = Node[]
    for e in node.edge
        if ischildof(node, e)
            push!(parents, getparent(e))
        end
    end
    return parents
end

@doc (@doc getparent) getparentedge
@inline function getparentedge(n::Node)
    for ee in n.edge
        if ee.isMajor && ischildof(n,ee)
            return ee
        end
    end
    error("node $(n.number) has no major parent")
end
@doc (@doc getparent) getparentedgeminor
@inline function getparentedgeminor(n::Node)
    for ee in n.edge
        if !ee.isMajor && n == ee.node[(ee.isChild1 ? 1 : 2)]
            return ee
        end
    end
    error("node $(n.number) has no minor parent")
end

"""
    getpartneredge(edge::Edge)
    getpartneredge(edge::Edge, node::Node)

Edge that is the hybrid partner of `edge`, meaning that is has the same child
`node` as `edge`. This child `node` is given as an argument in the second method.
Assumptions, not checked:

- no in-coming polytomy: a node has 0, 1 or 2 parents, no more
- when `node` is given, it is assumed to be the child of `edge`
  (the first method calls the second).

See also: [`getparent`](@ref), [`getchild`](@ref)
"""
@inline function getpartneredge(edge)
    node = getchild(edge)
    getpartneredge(edge, node)
end
@inline function getpartneredge(edge::Edge, node::Node)
    for e in node.edge
        if e.hybrid && e !== edge && node === getchild(e)
            return e
        end
    end
    error("did not find a partner for edge $(edge.number)")
end

"""
    edgerelation(e::Edge, node::Node, origin::Edge)

Return a symbol:

- `:origin` if `e` is equal to `origin`, and otherwise:
- `:parent` if `e` is a parent of `node`,
- `:child` if `e` is a child of `node`

using the `isChild1` attribute of edges.
Useful when `e` iterates over all edges adjacent to `node` and when
`origin` is one of the edges adjacent to `node`,
to known the order in which these edges come.

example:
```julia
labs = [edgerelation(e, u, uv) for e in u.edge] # assuming u is a node of edge uv
parentindex = findfirst(isequal(:parent), labs) # could be 'nothing' if no parent
childindices = findall( isequal(:child), labs)  # vector. could be empty
```
"""
function edgerelation(ee::Edge, n::Node, origin::Edge)
    (ee===origin ? :origin : (n===getchild(ee) ? :parent : :child))
end

# -------------- NODE -------------------------#

function setEdge!(node::Node,edge::Edge)
   push!(node.edge,edge);
   node.hasHybEdge = any(e -> e.hybrid, node.edge)
end

function getOtherNode(edge::Edge, node::Node)
  edge.node[1] === node ? edge.node[2] : edge.node[1]
end



# -------------- NETWORK ----------------------- #

function getIndex(node::Node, net::Network)
    i = 1;
    while(i<= size(net.node,1) && !isEqual(node,net.node[i]))
        i = i+1;
    end
    i <= size(net.node,1) || error("node $(node.number) not in network")
    return i
end

function getIndex(edge::Edge, net::Network)
    i = 1;
    while(i<= size(net.edge,1) && !isEqual(edge,net.edge[i]))
        i = i+1;
    end
    i <= size(net.edge,1) || error("edge $(edge.number) not in network")
    return i
end

function getIndex(edge::Edge, edges::Vector{Edge})
    i = 1;
    while(i<= size(edges,1) && !isEqual(edge,edges[i]))
        i = i+1;
    end
    i <= size(edges,1) || error("edge $(edge.number) not in array of edges")
    return i
end

# aux function to find the index of a node in a
# node array
function getIndex(name::Node, array::Array{Node,1})
    i = 1;
    while(i<= size(array,1) && !isequal(name,array[i]))
        i = i+1;
    end
    i <= size(array,1) || error("$(name.number) not in array")
    return i
end


function getIndexNode(number::Integer,net::Network)
    ind = findfirst(n -> n.number == number, net.node)
    if ind === nothing
        error("node number not in net.node")
    end
    return ind
end

function getIndexEdge(number::Integer,net::Network)
    ind = findfirst(x -> x.number == number, net.edge)
    if ind === nothing
        error("edge number not in net.edge")
    end
    return ind
end

# find the index of an edge in node.edge
function getIndexEdge(edge::Edge,node::Node)
    findfirst(e -> isequal(edge,e), node.edge)
end

# find the index of an edge with given number in node.edge
# bug found & fixed 2019-08-22. Unused function?
function getIndexEdge(number::Integer,node::Node)
    findfirst(e -> isequal(number,e.number), node.edge)
end

# find the index of a node in edge.node
function getIndexNode(edge::Edge,node::Node)
    size(edge.node,1) == 2 || @warn "this edge $(edge.number) has more or less than 2 nodes: $([n.number for n in edge.node])"
    if isequal(node,edge.node[1])
        return 1
    elseif isequal(node,edge.node[2])
        return 2
    else
        error("node not in edge.node")
    end
end

# function to find hybrid index in net.hybrid
function getIndexHybrid(node::Node, net::Network)
    node.hybrid || error("node $(node.number) is not hybrid so it cannot be in net.hybrid")
    i = 1;
    while(i<= size(net.hybrid,1) && !isEqual(node,net.hybrid[i]))
        i = i+1;
    end
    if i>size(net.hybrid,1) error("hybrid node not in network"); end
    return i
end

# function that given a hybrid node, it gives you the minor hybrid edge
# warning: assumes level-1 network: see getparentedgeminor for a general network
function getHybridEdge(node::Node)
    node.hybrid || error("node $(node.number) is not hybrid node, cannot get hybrid edges")
    a = nothing;
    for e in node.edge
        (e.hybrid && !e.isMajor) ? a = e : nothing; # assumes level-1: child of hybrid node must be a tree edge
    end
    isa(a,Nothing) ? error("hybrid node $(node.number) does not have minor hybrid edge, edges: $([e.number for e in node.edge])") : return a
end


# function that given two nodes, it gives you the edge that connects them
# returns error if they are not connected by an edge
function getConnectingEdge(node1::Node,node2::Node)
    found = false;
    i = 1;
    while(i<= size(node1.edge,1) && !found)
        if(isequal(getOtherNode(node1.edge[i],node1),node2))
            found = true;
        end
        i = i+1;
    end
    if(found)
        return node1.edge[i-1]
    else
        error("nodes not connected")
    end
end

"""
    isconnected(node1::Node, node2::Node)

Check if two nodes are connected by an edge. Return true if connected, false
if not connected.
"""
function isconnected(node1, node2)
    for e in node1.edge
        if e in node2.edge
            return true
        end
    end
    return false
end

# function to check in an edge is in an array by comparing
# edge numbers (could use isEqual for adding comparisons of gammaz and inCycle)
# needed for updateHasEdge
function isEdgeNumIn(edge::Edge,array::Array{Edge,1})
    enum = edge.number
    return any(e -> e.number == enum, array)
end

# function to check in a leaf is in an array by comparing
# the numbers (uses isEqual)
# needed for updateHasEdge
function isNodeNumIn(node::Node,array::Array{Node,1})
    return all((e->!isEqual(node,e)), array) ? false : true
end

# function to push a Node in net.node and
# update numNodes and numTaxa
function pushNode!(net::Network, n::Node)
    push!(net.node,n);
    net.numNodes += 1;
    if(n.leaf)
        net.numTaxa += 1
        push!(net.leaf,n);
    end
    if(n.hybrid)
        pushHybrid!(net,n)
    end
end

# function to push an Edge in net.edge and
# update numEdges
function pushEdge!(net::Network, e::Edge)
    push!(net.edge,e);
    net.numEdges += 1;
end


# function to push a hybrid Node in net.hybrid and
# update numHybrids
function pushHybrid!(net::Network, n::Node)
    if(n.hybrid)
        push!(net.hybrid,n);
        net.numHybrids += 1;
    else
        error("node $(n.number) is not hybrid, so cannot be pushed in net.hybrid")
    end
end

"""
    deleteNode!(net::HybridNetwork, n::Node)
    deleteNode!(net::QuartetNetwork, n::Node)

Delete node `n` from a network, i.e. removes it from
net.node, and from net.hybrid or net.leaf as appropriate.
Update attributes `numNodes`, `numTaxa`, `numHybrids`.
Warning: `net.names` is *not* updated, and this is a feature (not a bug)
for networks of type QuartetNetwork.

Warning: if the root is deleted, the new root is arbitrarily set to the
first node in the list. This is intentional to save time because this function
is used frequently in snaq!, which handles semi-directed (unrooted) networks.
"""

function deleteNode!(net::HybridNetwork, n::Node;thorough=false::Bool)
    index = 0
    try
        ##Note from Josh
            ##With the things I'm doing some nodes on a tree may have the same attributes causing getIndex() to return the wrong node
            ##Seems like using == to compare the nodes gets me the right thing. So I am manually getting the index here using ==. 
            ##Perhaps a more appropriate place to have a flag for the type of comparison is in the getIndex function itself. TODO 
        if(!thorough)
            index = getIndex(n,net);
        else ##Get the index manually with ==
            index = findfirst(x -> x===n, net.node)
        end     
    catch
        error("Node $(n.number) not in network");
    end
	index !== nothing || error("Node $(n.number) not in network"); ##Might be redundant
    # println("deleting node $(n.number) from net, index $(index).")
    deleteat!(net.node,index);
    net.numNodes -= 1;
    if net.root == index  # do not check containRoot to save time in snaq!
        net.root = 1      # arbitrary
    elseif net.root > index
        net.root -= 1
    end
    if n.hybrid
       removeHybrid!(net,n)
    end
    if n.leaf
        removeLeaf!(net,n)
    end
end

function deleteNode!(net::QuartetNetwork, n::Node)
    index = findfirst(no -> no.number == n.number, net.node)
    # isEqual (from above) checks for more than node number
    index !== nothing || error("Node $(n.number) not in quartet network");
    deleteat!(net.node,index);
    net.numNodes -= 1
    if n.hybrid
       removeHybrid!(net,n)
    end
    if n.leaf
        index = findfirst(no -> no === n, net.leaf)
        index !== nothing || error("node $(n.number) not net.leaf")
        deleteat!(net.leaf,index)
        net.numTaxa -= 1
    end
end

"""
    deleteEdge!(net::HybridNetwork,  e::Edge; part=true)
    deleteEdge!(net::QuartetNetwork, e::Edge)

Delete edge `e` from `net.edge` and update `net.numEdges`.
If `part` is true, update the network's partition field.
"""
function deleteEdge!(net::HybridNetwork, e::Edge; part=true::Bool)
    if part
        if e.inCycle == -1 && !e.hybrid && !isempty(net.partition) && !isTree(net)
            ind = whichPartition(net,e)
            indE = getIndex(e,net.partition[ind].edges)
            deleteat!(net.partition[ind].edges,indE)
        end
    end
    i = findfirst(x -> x===e, net.edge)
    i !== nothing || error("edge $(e.number) not in network: can't delete");
    deleteat!(net.edge, i);
    net.numEdges -= 1;
end

# function to delete an Edge in net.edge and
# update numEdges from a QuartetNetwork
function deleteEdge!(net::QuartetNetwork, e::Edge)
    index = findfirst(x -> x.number == e.number, net.edge)
    # isEqual (from above) checks for more than edge number
    index !== nothing || error("edge not in quartet network");
    deleteat!(net.edge,index);
    net.numEdges -= 1;
end


"""
    removeHybrid!(net::Network, n::Node)

Delete a hybrid node `n` from `net.hybrid`, and update `net.numHybrid`.
The actual node `n` is not deleted. It is kept in the full list `net.node`.
"""
function removeHybrid!(net::Network, n::Node)
    n.hybrid || error("cannot delete node $(n.number) from net.hybrid because it is not hybrid")
    i = findfirst(x -> x===n, net.hybrid)
    i !== nothing || error("hybrid node $(n.number) not in the network's list of hybrids");
    deleteat!(net.hybrid, i);
    net.numHybrids -= 1;
end

# function to delete a leaf node in net.leaf
# and update numTaxa
function removeLeaf!(net::Network,n::Node)
    n.leaf || error("cannot delete node $(n.number) from net.leaf because it is not leaf")
    index = findfirst(no -> no === n, net.leaf)
    index !== nothing || error("leaf node $(n.number) not in network")
    deleteat!(net.leaf,index)
    net.numTaxa -= 1
end

# function to delete an internal node with only 2 edges
function deleteIntNode!(net::Network, n::Node)
    size(n.edge,1) == 2 || error("node $(n.number) does not have only two edges")
#    isEqual(n,net.node[net.root]) && println("deleting the root $(n.number) because it has only two edges attached")
    index = n.edge[1].number < n.edge[2].number ? 1 : 2;
    edge1 = n.edge[index];
    edge2 = n.edge[index==1 ? 2 : 1];
    if(!edge1.hybrid && !edge2.hybrid)
        node1 = getOtherNode(edge1,n);
        node2 = getOtherNode(edge2,n);
        removeEdge!(node2,edge2);
        removeNode!(n,edge1);
        setEdge!(node2,edge1);
        setNode!(edge1,node2);
        deleteNode!(net,n);
        deleteEdge!(net,edge2);
    else
        @warn "the two edges $([edge1.number,edge2.number]) attached to node $(n.number) must be tree edges to delete node"
        if(edge1.hybrid)
            hybedge = edge1
            otheredge = edge2
        elseif(edge2.hybrid)
            hybedge = edge2
            otheredge = edge1
        end
        othernode = getOtherNode(otheredge,n)
        removeNode!(n,hybedge)
        removeEdge!(othernode,otheredge)
        setEdge!(othernode,hybedge)
        setNode!(hybedge,othernode)
        deleteNode!(net,n)
        deleteEdge!(net,otheredge)
    end
end


# search the hybrid node(s) in network: returns the hybrid node(s)
# in an array
# throws error if no hybrid in network
function searchHybridNode(net::Network)
    a = [n for n in net.node if n.hybrid]
    suma = length(a)
    suma != 0 || error("network has no hybrid node")
    return a
end

# search and returns the hybrid edges in network
# throws error if no hybrid in network
function searchHybridEdge(net::Network)
    a = [n for n in net.edge if n.hybrid]
    suma = length(a)
    suma != 0 || error("network has no hybrid edge")
    return a
end

"""
    printEdges(net)
    printEdges(io::IO, net)

Print information on the edges of a `HybridNetwork` or `QuartetNetwork` object
`net`: edge number, numbers of nodes attached to it, edge length, whether it's
a hybrid edge, its γ inheritance value, whether it's a major edge,
if it could contain the root (this field is not always updated, though)
and attributes pertaining to level-1 networks used in SNaQ:
in which cycle it is contained (-1 if no cycle), and if the edge length
is identifiable (based on quartet concordance factors).
"""
printEdges(x) = printEdges(stdout::IO, x)
function printEdges(io::IO, net::HybridNetwork)
    if net.numBad > 0
        println(io, "net has $(net.numBad) bad diamond I. Some γ and edge lengths t are not identifiable, although their γ * (1-exp(-t)) are.")
    end
    miss = ""
    println(io, "edge parent child  length  hybrid isMajor gamma   containRoot inCycle istIdentitiable")
    for e in net.edge
        @printf(io, "%-4d %-6d %-6d ", e.number, getparent(e).number, getchild(e).number)
        if e.length==-1.0 @printf(io, "%-7s ", miss); else @printf(io, "%-7.3f ", e.length); end
        @printf(io, "%-6s %-7s ", e.hybrid, e.isMajor)
        if e.gamma==-1.0  @printf(io, "%-7s ", miss); else @printf(io, "%-7.4g ", e.gamma); end
        @printf(io, "%-11s %-7d %-5s\n", e.containRoot, e.inCycle, e.istIdentifiable)
    end
end

function printEdges(io::IO, net::QuartetNetwork)
    println(io, "edge parent child  length  hybrid isMajor gamma   containRoot inCycle istIdentitiable")
    for e in net.edge
        @printf(io, "%-4d %-6d %-6d ", e.number, getparent(e).number, getchild(e).number)
        @printf(io, "%-7.3f %-6s %-7s ", e.length, e.hybrid, e.isMajor)
        @printf(io, "%-7.4g %-11s %-7d %-5s\n", e.gamma, e.containRoot, e.inCycle, e.istIdentifiable)
    end
end

# print for every node, inCycle and edges
"""
    printNodes(net)
    printNodes(io, net)

Print information on the nodes of a `HybridNetwork` net: node number,
whether it's a leaf, whether it's a hybrid node, whether it's connected to one
or more hybrid edges, it's name (label),
the cycle in which it is belong (-1 if no cycle; makes sense for level-1 networks),
and the list of edges attached to it, by their numbers.
"""
printNodes(x) = printNodes(stdout::IO, x)
function printNodes(io::IO, net::Network)
    namepad = max(4, maximum(length.([n.name for n in net.node])))
    println(io, "node leaf  hybrid hasHybEdge ", rpad("name", namepad), " inCycle edges'numbers")
    for n in net.node
        @printf(io, "%-4d %-5s %-6s %-10s ", n.number, n.leaf, n.hybrid, n.hasHybEdge)
        print(io, rpad(n.name,namepad))
        @printf(io, " %-7d", n.inCycle)
        for e in n.edge
            @printf(io, " %-4d", e.number)
        end
        print(io, "\n")
    end
end

"""
    hybridEdges(node::Node)

Return the 3 edges attached to `node` in a specific order [e1,e2,e3].
**Warning**: assume a level-1 network with node field `hasHybEdge`
and edge field `inCycle` up-to-date.

If `node` is a hybrid node:

- e1 is the major hybrid parent edge of `node`
- e2 is the minor hybrid parent edge
- e3 is the tree edge, child of `node`.

If `node` is a tree node parent of one child edge:

- e1 is the hybrid edge, child of `node`
- e2 is the tree edge that belongs to the cycle created by e1
- e3 is the other tree edge attached to `node` (not in a cycle)

Otherwise:

- e3 is an external edge from `node` to a leaf, if one exists.
"""
function hybridEdges(node::Node)
    size(node.edge,1) == 3 || error("node $(node.number) has $(size(node.edge,1)) edges instead of 3");
    if(node.hybrid)
        hybmajor = nothing;
        hybminor = nothing;
        tree = nothing;
        for e in node.edge
            (e.hybrid && e.isMajor) ? hybmajor = e : nothing
            (e.hybrid && !e.isMajor) ? hybminor = e : nothing
            !e.hybrid ? tree = e : nothing
        end
        return hybmajor, hybminor, tree
    elseif(node.hasHybEdge)
        hybrid = nothing;
        treecycle = nothing;
        tree = nothing;
        for e in node.edge
            (e.hybrid) ? hybrid = e : nothing
            (!e.hybrid && e.inCycle != -1) ? treecycle = e : nothing
            (!e.hybrid && e.inCycle == -1) ? tree = e : nothing
        end
        return hybrid, treecycle, tree
    else
        #@warn "node $(node.number) is not hybrid $(node.hybrid) nor tree with hybrid edges (hasHybEdge) $(node.hasHybEdge), return the node.edge in order, unless a leaf is attached, then the edge attached to leaf is last";
        edge1 = nothing
        edge2 = nothing
        edge3 = nothing
        leaffound = false
        ind = 1
        for i in 1:3
            if(getOtherNode(node.edge[i],node).leaf)
                leaffound = true
                edge3 = node.edge[i]
                ind = i
                break
            end
        end
        if(leaffound)
            if(ind == 1)
                return node.edge[2], node.edge[3], edge3
            elseif(ind == 2)
                return node.edge[1], node.edge[3], edge3
            elseif(ind == 3)
                return node.edge[1], node.edge[2], edge3
            end
        else
            return node.edge[1], node.edge[2], node.edge[3]
        end
    end
end

"""
    hybridEdges(node::Node, e::Edge)

Return the 2 edges connected to `node` other than `e`,
in the same order as `node.edge`,
except that `e` absent from the list.

Despite what the name suggest, `node` need not be a hybrid node!
`node` is assumed to have 3 edges, though.
"""
function hybridEdges(node::Node, edge::Edge)
    size(node.edge,1) == 3 || error("node $(node.number) has $(size(node.edge,1)) edges instead of 3")
    edge1 = nothing
    edge2 = nothing
    for e in node.edge
        if(!isequal(e,edge))
            isa(edge1,Nothing) ? edge1 = e : edge2 = e
        end
    end
    return edge1,edge2
end


# function to remove an edge from a node
# warning: deletion is final, you can only
#          have edge back by pushing it again
# warning: if the edge removed is hybrid and node is tree,
#          node.hasHybEdge is set to false
#          assuming any tree node can only have one
#          one hybrid edge
function removeEdge!(node::Node, edg::Edge)
    index = findfirst(x -> x === edg, node.edge)
    index !== nothing || error("edge $(edg.number) not in node $(node.number)")
    deleteat!(node.edge,index)
    node.hasHybEdge = any(e -> e.hybrid, node.edge)
end

# function to remove a node from a edge
# warning: deletion is final, you can only
#          have node back by pushing it again
# warning: only removes node from edge, edge might still
#          be in node.edge
function removeNode!(nod::Node, edge::Edge)
    index = findfirst(x -> x === nod, edge.node)
    index !== nothing || error("node $(nod.number) not in edge")
    deleteat!(edge.node,index);
end


# ----------------------------------------------------------------------------------------

# setLength
# warning: allows to change edge length for istIdentifiable=false
#          but issues a warning
# negative=true means it allows negative branch lengths (useful in qnet typeHyb=4)
function setLength!(edge::Edge, new_length::Number, negative::Bool)
    (negative || new_length >= 0) || error("length has to be nonnegative: $(new_length), cannot set to edge $(edge.number)")
    new_length >= -0.4054651081081644 || error("length can be negative, but not too negative (greater than -log(1.5)) or majorCF<0: new length is $(new_length)")
    #println("setting length $(new_length) to edge $(edge.number)")
    if(new_length > 10.0)
        new_length = 10.0;
    end
    edge.length = new_length;
    edge.y = exp(-new_length);
    edge.z = 1.0 - edge.y;
    #edge.istIdentifiable || @warn "set edge length for edge $(edge.number) that is not identifiable"
    return nothing
end

"""
    setLength!(edge, newlength)`

Set the length of `edge`, and set `edge.y` and `edge.z` accordingly.
Warning: specific to SNaQ. Use [`setlengths!`](@ref) or [`setBranchLength!`](@ref)
for more general tools.

- The new length is censored to 10: if the new length is above 10,
  the edge's length will be set to 10. Lengths are interpreted in coalescent
  units, and 10 is close to infinity: near perfect gene tree concordance.
  10 is used as an upper limit to coalescent units that can be reliably estimated.
- The new length is allowed to be negative, but must be greater than -log(1.5),
  to ensure that the major quartet concordance factor (1 - 2/3 exp(-length)) is >= 0.
"""
setLength!(edge::Edge, new_length::Number) = setLength!(edge, new_length, false)


"""
    setBranchLength!(Edge, newlength)

Set the length of an Edge object. The new length needs to be non-negative,
or -1.0 to be interpreted as missing. `edge.y` and `edge.z` are updated
accordingly.
"""
function setBranchLength!(edge::Edge, new_length::Number)
    (new_length >= 0 || new_length == -1.0) || error("length $(new_length) has to be nonnegative or -1.0 (for missing).")
    edge.length = new_length;
    edge.y = exp(-new_length);
    edge.z = 1.0 - edge.y;
end


"""
    setGamma!(Edge, new γ)
    setGamma!(Edge, new γ, change_other=true::Bool)

Set inheritance probability γ for an edge, which must be a hybrid edge.
The new γ needs to be in [0,1]. The γ of the "partner" hybrid edge is changed
accordingly, to 1-γ. The field `isMajor` is also changed accordingly.
If the new γ is approximately 0.5, `Edge` is set to the major parent,
its partner is set to the minor parent.

If `net` is a HybridNetwork object, `printEdges(net)` will show the list of edges
and their γ's. The γ of the third hybrid edge (say) can be changed to 0.2 with
`setGamma!(net.edge[3],0.2)`.
This will automatically set γ of the partner hybrid edge to 0.8.

The last argument is true by default. If false: the partner edge is not updated.
This is useful if the new γ is 0.5, and the partner's γ is already 0.5,
in which case the `isMajor` attributes can remain unchanged.
"""
setGamma!(edge::Edge, new_gamma::Float64) = setGamma!(edge, new_gamma, true)

# warning in the bad diamond/triangle cases because gamma is not identifiable
# changeOther = true, looks for the other hybrid edge and changes gamma too

function setGamma!(edge::Edge, new_gamma::Float64, changeOther::Bool)
    new_gamma >= 0.0 || error("gamma has to be positive: $(new_gamma)")
    new_gamma <= 1.0 || error("gamma has to be less than 1: $(new_gamma)")
    edge.hybrid || error("cannot change gamma in a tree edge");
    node = getchild(edge) # child of hybrid edge
    node.hybrid || @warn "hybrid edge $(edge.number) not pointing at hybrid node"
    # @debug (node.isBadDiamondI ? "bad diamond situation: gamma not identifiable" : "")
    partner = Edge[] # list of other hybrid parents of node, other than edge
    for e in node.edge
        if e.hybrid && e != edge && node == getchild(e)
            push!(partner, e)
        end
    end
    length(partner) == 1 ||
      error("strange hybrid node $(node.number) with $(length(partner)+1) hybrid parents")
    e2 = partner[1]
    onehalf = isapprox(new_gamma,0.5)
    if onehalf new_gamma=0.5; end
    new_ismajor = new_gamma >= 0.5
    edge.gamma = new_gamma
    if changeOther
        edge.isMajor = new_ismajor
        e2.gamma = 1.0 - new_gamma
        e2.isMajor = !new_ismajor
    else
        if onehalf # who is major is arbitrary: so we pick what's consistent with the partner
            edge.isMajor = !e2.isMajor
        else
            edge.isMajor = new_ismajor
        end
    end
    return nothing
end

@inline function setmultiplegammas!(edges::Vector{Edge}, gammas::Vector{Float64})
    for (e,g) in zip(edges, gammas)
        setGamma!(e, g)
    end
end

"""
    setGammaBLfromGammaz!(node, network)

Update the γ values of the two sister hybrid edges in a bad diamond I, given the `gammaz` values
of their parent nodes, and update the branch lengths t1 and t2 of their parent edges
(those across from the hybrid nodes), in such a way that t1=t2 and that these branch lengths
and γ values are consistent with the `gammaz` values in the network.

Similar to the first section of [`undoGammaz!`](@ref),
but does not update anything else than γ and t's.
Unlike `undoGammaz!`, no error if non-hybrid `node` or not at bad diamond I.
"""
function setGammaBLfromGammaz!(node::Node, net::HybridNetwork)
    if !node.isBadDiamondI || !node.hybrid
        return nothing
    end
    edge_maj, edge_min, tree_edge2 = hybridEdges(node);
    other_maj = getOtherNode(edge_maj,node);
    other_min = getOtherNode(edge_min,node);
    edgebla,tree_edge_incycle1,tree_edge = hybridEdges(other_min);
    edgebla,tree_edge_incycle2,tree_edge = hybridEdges(other_maj);
    if(approxEq(other_maj.gammaz,0.0) && approxEq(other_min.gammaz,0.0))
        edge_maj.gamma = 1.0 # γ and t could be anything if both gammaz are 0
        edge_min.gamma = 0.0 # will set t's to 0 and minor γ to 0.
        newt = 0.0
    else
        ((approxEq(other_min.gammaz,0.0) || other_min.gammaz >= 0.0) &&
         (approxEq(other_maj.gammaz,0.0) || other_maj.gammaz >= 0.0)    ) ||
            error("bad diamond I in node $(node.number) but missing (or <0) gammaz")
        ztotal = other_maj.gammaz + other_min.gammaz
        edge_maj.gamma = other_maj.gammaz / ztotal
        edge_min.gamma = other_min.gammaz / ztotal
        newt = -log(1-ztotal)
    end
    setLength!(tree_edge_incycle1,newt)
    setLength!(tree_edge_incycle2,newt)
end


function numTreeEdges(net::HybridNetwork)
    2*net.numTaxa - 3 + net.numHybrids
end

function numIntTreeEdges(net::HybridNetwork)
    2*net.numTaxa - 3 + net.numHybrids - net.numTaxa
end


# function to get the partition where an edge is
# returns the index of the partition, or error if not found
# better to return the index than the partition itself, because we need the index
# to use splice and delete it from net.partition later on
# cycle: is the number to look for partition on that cycle only
function whichPartition(net::HybridNetwork,edge::Edge,cycle::Integer)
    !edge.hybrid || error("edge $(edge.number) is hybrid so it cannot be in any partition")
    edge.inCycle == -1 || error("edge $(edge.number) is in cycle $(edge.inCycle) so it cannot be in any partition")
    @debug "search partition for edge $(edge.number) in cycle $(cycle)"
    in(edge,net.edge) || error("edge $(edge.number) is not in net.edge")
    for i in 1:length(net.partition)
        @debug "looking for edge $(edge.number) in partition $(i): $([e.number for e in net.partition[i].edges])"
        if(in(cycle,net.partition[i].cycle))
            @debug "looking for edge $(edge.number) in partition $(i), with cycle $(cycle): $([e.number for e in net.partition[i].edges])"
            if in(edge,net.partition[i].edges)
                @debug "partition for edge $(edge.number) is $([e.number for e in net.partition[i].edges])"
                return i
            end
        end
    end
    @debug begin; printPartitions(net); "" end
    error("edge $(edge.number) is not hybrid, nor part of any cycle, and it is not in any partition")
end

# function to get the partition where an edge is
# returns the index of the partition, or error if not found
# better to return the index than the partition itself, because we need the index
# to use splice and delete it from net.partition later on
function whichPartition(net::HybridNetwork,edge::Edge)
    !edge.hybrid || error("edge $(edge.number) is hybrid so it cannot be in any partition")
    edge.inCycle == -1 || error("edge $(edge.number) is in cycle $(edge.inCycle) so it cannot be in any partition")
    @debug "search partition for edge $(edge.number) without knowing its cycle"
    in(edge,net.edge) || error("edge $(edge.number) is not in net.edge")
    for i in 1:length(net.partition)
        @debug "looking for edge $(edge.number) in partition $(i): $([e.number for e in net.partition[i].edges])"
        if(in(edge,net.partition[i].edges))
            @debug "partition for edge $(edge.number) is $([e.number for e in net.partition[i].edges])"
            return i
        end
    end
    @debug begin printPartitions(net); "printed partitions" end
    error("edge $(edge.number) is not hybrid, nor part of any cycle, and it is not in any partition")
end

# function that will print the partition of net
function printPartitions(net::HybridNetwork)
    println("partition.cycle\t partition.edges")
    for p in net.partition
        println("$(p.cycle)\t\t $([e.number for e in p.edges])")
    end
end

# function to find if a given partition is in net.partition
function isPartitionInNet(net::HybridNetwork,desc::Vector{Edge},cycle::Vector{Int})
    for p in net.partition
        if(sort(cycle) == sort(p.cycle))
            if(sort([e.number for e in desc]) == sort([e.number for e in p.edges]))
                return true
            end
        end
    end
    return false
end

# function to check that everything matches in a network
# in particular, cycles, partitions and containRoot
# fixit: need to add check on identification of bad diamonds, triangles
# and correct computation of gammaz
# light=true: it will not collapse with nodes with 2 edges, will return a flag of true
# returns true if found egde with BL -1.0 (only when light=true, ow error)
# added checkPartition for undirectedOtherNetworks that do not need correct hybrid node number
function checkNet(net::HybridNetwork, light::Bool; checkPartition=true::Bool)
    @debug "checking net"
    net.numHybrids == length(net.hybrid) || error("discrepant number on net.numHybrids (net.numHybrids) and net.hybrid length $(length(net.hybrid))")
    net.numTaxa == length(net.leaf) || error("discrepant number on net.numTaxa (net.numTaxa) and net.leaf length $(length(net.leaf))")
    net.numNodes == length(net.node) || error("discrepant number on net.numNodes (net.numNodes) and net.node length $(length(net.node))")
    net.numEdges == length(net.edge) || error("discrepant number on net.numEdges (net.numEdges) and net.edge length $(length(net.edge))")
    if(isTree(net))
        all(x->x.containRoot,net.edge) || error("net is a tree, but not all edges can contain root")
        all(x->x.isMajor,net.edge) || error("net is a tree, but not all edges are major")
        all(x->!(x.hybrid),net.edge) || error("net is a tree, but not all edges are tree")
        all(x->!(x.hybrid),net.node) || error("net is a tree, but not all nodes are tree")
        all(x->!(x.hasHybEdge),net.node) || error("net is a tree, but not all nodes hasHybEdge=false")
        all(x->(x.gamma == 1.0 ? true : false),net.edge) || error("net is a tree, but not all edges have gamma 1.0")
    end
    for h in net.hybrid
        if(isBadTriangle(h))
            @debug "hybrid $(h.number) is very bad triangle"
            net.hasVeryBadTriangle || error("hybrid node $(h.number) is very bad triangle, but net.hasVeryBadTriangle is $(net.hasVeryBadTriangle)")
            h.isVeryBadTriangle || h.isExtBadTriangle || error("hybrid node $(h.number) is very bad triangle but it does not know it")
        end
        nocycle,edges,nodes = identifyInCycle(net,h)
        for e in edges
            e.inCycle == h.number || error("edge $(e.number) is in cycle of hybrid node $(h.number) but its inCycle attribute is $(e.inCycle)")
            if(e.length == -1.0)
                if(light)
                    return true
                else
                    error("found edge with BL -1.0")
                end
            end
            if(e.hybrid)
                !e.containRoot || error("hybrid edge $(e.number) should not contain root") # fixit: disagree
                o = getOtherNode(e,h)
                o.hasHybEdge || error("found node $(o.number) attached to hybrid edge but hasHybEdge=$(o.hasHybEdge)")
            end
        end
        for n in nodes
            n.inCycle == h.number || error("node $(n.number) is in cycle of hybrid node $(h.number) but its inCycle attribute is $(n.inCycle)")
            e1,e2,e3 = hybridEdges(n)
            i = 0
            for e in [e1,e2,e3]
                if(isa(e,Nothing) && h.k != 2)
                    error("edge found that is Nothing, and hybrid node $(h.number) k is $(h.k). edge as nothing can only happen when k=2")
                elseif(!isa(e,Nothing))
                    if(e.inCycle == -1)
                        i += 1
                        desc = [e]
                        cycleNum = [h.number]
                        getDescendants!(getOtherNode(e,n),e,desc,cycleNum)
                        if(checkPartition && !isPartitionInNet(net,desc,cycleNum))
                            printPartitions(net)
                            error("partition with cycle $(cycleNum) and edges $([e.number for e in desc]) not found in net.partition")
                        end
                    end
                end
            end
            i == 1 || error("strange node $(n.number) incycle $(h.number) but with $(i) edges not in cycle, should be only one")
            edgesRoot = identifyContainRoot(net,h)
            for edge in edgesRoot
                if edge.containRoot
                    @debug begin printEverything(net); "printed everything" end
                    error("edge $(edge.number) should not contain root")
                end
            end
        end
    end
    for n in net.node
        if(n.leaf)
            length(n.edge) == 1 || error("leaf $(n.number) with $(length(n.edge)) edges instead of 1")
        else
            if(light)
                if(length(n.edge) != 3)
                    @debug "warning: node $(n.number) with $(length(n.edge)) edges instead of 3"
                    return true
                end
            else
                length(n.edge) == 3 || error("node $(n.number) with $(length(n.edge)) edges instead of 3")
            end
        end
    end
    for e in net.edge
        if(e.length == -1.0)
            if(light)
                return true
            else
                error("edge found with BL -1.0")
            end
        end
    end
    @debug "no errors in checking net"
    return false
end

checkNet(net::HybridNetwork) = checkNet(net, false)

# function to print everything for a given net
# this is used a lot inside snaq to debug, so need to use level1 attributes
# and not change the network: with writeTopologyLevel1
function printEverything(net::HybridNetwork)
    printEdges(net)
    printNodes(net)
    printPartitions(net)
    println("$(writeTopologyLevel1(net))")
end

# function to check if a node is very or ext bad triangle
function isBadTriangle(node::Node)
    node.hybrid || error("cannot check if node $(node.number) is very bad triangle because it is not hybrid")
    if(node.k == 3)
        edgemaj, edgemin, treeedge = hybridEdges(node)
        othermaj = getOtherNode(edgemaj,node)
        othermin = getOtherNode(edgemin,node)
        treenode = getOtherNode(treeedge,node)
        edges1 = hybridEdges(othermaj)
        o1 = getOtherNode(edges1[3],othermaj)
        edges2 = hybridEdges(othermin)
        o2 = getOtherNode(edges2[3],othermin)
        leaves = sum([n.leaf ? 1 : 0 for n in [treenode,o1,o2]])
        if(leaves == 1 || leaves == 2)
            return true
        else
            return false
        end
    else
        return false
    end
end


# function to check if a partition is already in net.partition
# used in updatePartition
function isPartitionInNet(net::HybridNetwork,partition::Partition)
    if(isempty(net.partition))
        return false
    end
    for p in net.partition
        cycle = isempty(setdiff(p.cycle,partition.cycle)) && isempty(setdiff(partition.cycle,p.cycle))
        edges = isempty(setdiff([n.number for n in p.edges],[n.number for n in partition.edges])) && isempty(setdiff([n.number for n in partition.edges],[n.number for n in p.edges]))
        if(cycle && edges)
            return true
        end
    end
    return false
end


# function to switch a hybrid node in a network to another node in the cycle
function switchHybridNode!(net::HybridNetwork, hybrid::Node, newHybrid::Node)
    hybrid.hybrid || error("node $(hybrid.number) has to be hybrid to switch to a different hybrid")
    newHybrid.inCycle == hybrid.number || error("new hybrid needs to be in the cycle of old hybrid: $(hybrid.number)")
    !newHybrid.hybrid || error("strange hybrid node $(newHybrid.number) in cycle of another hybrid $(hybrid.number)")
    newHybrid.hybrid = true
    newHybrid.hasHybEdge = true
    newHybrid.name = hybrid.name
    pushHybrid!(net,newHybrid)
    makeNodeTree!(net,hybrid)
end

"""
    assignhybridnames!(net)

Assign names to hybrid nodes in the network `net`.
Hybrid nodes with an empty `name` field ("") are modified with a name that
does not conflict with other hybrid names in the network. The preferred name
is "H3" if the node number is 3 or -3, but an index other than 3 would be used
if "H3" were the name of another node already.

If two hybrid nodes have non-empty and equal names, the name of one of them is changed and
re-assigned as described above (with a warning).
"""
function assignhybridnames!(net::HybridNetwork)
    rx = r"^H(\d+)$"
    # prep: collect indices 'i' of any tree nodes named like Hi
    trenum = Int[]  # indices 'i' in tree node name, in case some are named Hi
    for n in net.node
        !n.hybrid || continue # do nothing if node n is hybrid
        m = match(rx, n.name)
        m === nothing || push!(trenum, parse(Int, m[1]))
    end
    # first: go through *all* existing non-empty names
    hybnum = Int[]  # indices 'i' in hybrid names: Hi
    for ih in 1:length(net.hybrid)
        hnode = net.hybrid[ih]
        lab = hnode.name
        lab != "" || continue # do nothing if label is missing
        jh = findfirst(isequal(lab), [net.hybrid[j].name for j in 1:ih-1])
        if jh !== nothing # set repeated names to ""
            @warn "hybrid nodes $(hnode.number) and $(net.hybrid[jh].number) have the same label: $lab. Will change the name of the former."
            hnode.name = ""
        else # fill in list of existing indices "i" in Hi
            m = match(rx, lab)
            m !== nothing || continue # skip the rest if name is not of the form Hi
            ind = parse(Int, m[1])
            if ind in trenum
                @warn "hybrid node $(hnode.number) had same label as a tree node: H$ind. Will change hybrid name."
                hnode.name = ""
            else
                push!(hybnum, ind)
            end
        end
    end
    # second: assign empty names to "Hi" for some i
    hnext = 1
    for ih in 1:length(net.hybrid)
        net.hybrid[ih].name == "" || continue # do nothing if non-empty label
        hnum = abs(net.hybrid[ih].number)
        while hnum in hybnum || hnum in trenum
            hnum = hnext  # not efficient, but rare
            hnext += 1    # and okay on small networks
        end
        push!(hybnum, hnum)
        net.hybrid[ih].name = "H$hnum"
    end
end

"""
    sorttaxa!(DataFrame, columns)

Reorder the 4 taxa and reorders the observed concordance factors accordingly, on each row of
the data frame. If `columns` is ommitted, taxon names are assumed to be in columns 1-4 and
CFs are assumed to be in columns 5-6 with quartets in this order: `12_34`, `13_24`, `14_23`.
Does **not** reorder credibility interval values, if present.

    sorttaxa!(DataCF)
    sorttaxa!(Quartet, permutation_tax, permutation_cf)

Reorder the 4 taxa in each element of the DataCF `quartet`. For a given Quartet,
reorder the 4 taxa in its fields `taxon` and `qnet.quartetTaxon` (if non-empty)
and reorder the 3 concordance values accordingly, in `obsCF` and `qnet.expCF`.

`permutation_tax` and `permutation_cf` should be vectors of short integers (Int8) of length 4 and 3
respectively, whose memory allocation gets reused. Their length is *not checked*.

`qnet.names` is unchanged: the order of taxon names here relates to the order of nodes in the network
(???)
"""
function sorttaxa!(dat::DataCF)
    ptax = Array{Int8}(undef, 4) # to hold the sort permutations
    pCF  = Array{Int8}(undef, 3)
    for q in dat.quartet
        sorttaxa!(q, ptax, pCF)
    end
end

function sorttaxa!(df::DataFrame, co=Int[]::Vector{Int})
    if length(co)==0
        co = collect(1:7)
    end
    length(co) > 6 || error("column vector must be of length 7 or more")
    ptax = Array{Int8}(undef, 4)
    pCF  = Array{Int8}(undef, 3)
    taxnam = Array{eltype(df[!,co[1]])}(undef, 4)
    for i in 1:size(df,1)
        for j=1:4 taxnam[j] = df[i,co[j]]; end
        sortperm!(ptax, taxnam)
        sorttaxaCFperm!(pCF, ptax) # update permutation pCF according to taxon permutation
        df[i,co[1]], df[i,co[2]], df[i,co[3]], df[i,co[4]] = taxnam[ptax[1]], taxnam[ptax[2]], taxnam[ptax[3]], taxnam[ptax[4]]
        df[i,co[5]], df[i,co[6]], df[i,co[7]] = df[i,co[pCF[1]+4]], df[i,co[pCF[2]+4]], df[i,co[pCF[3]+4]]
    end
    return df
end

function sorttaxa!(qua::Quartet, ptax::Vector{Int8}, pCF::Vector{Int8})
    qt = qua.taxon
    if length(qt)==4
        sortperm!(ptax, qt)
        sorttaxaCFperm!(pCF, ptax) # update permutation pCF accordingly
        qt[1], qt[2], qt[3], qt[4] = qt[ptax[1]], qt[ptax[2]], qt[ptax[3]], qt[ptax[4]]
        qua.obsCF[1], qua.obsCF[2], qua.obsCF[3] = qua.obsCF[pCF[1]], qua.obsCF[pCF[2]], qua.obsCF[pCF[3]]
        # do *NOT* modify qua.qnet.quartetTaxon: it points to the same array as qua.taxon
        eCF = qua.qnet.expCF
        if length(eCF)==3
            eCF[1], eCF[2], eCF[3] = eCF[pCF[1]], eCF[pCF[2]], eCF[pCF[3]]
        end
    elseif length(qt)!=0
        error("Quartet with $(length(qt)) taxa")
    end
    return qua
end

# find permutation pCF of the 3 CF values: 12_34, 13_24, 14_23. 3!=6 possible permutations
# ptax = one of 4!=24 possible permutations on the 4 taxon names
# kernel: pCF = identity if ptax = 1234, 2143, 3412 or 4321
# very long code, but to minimize equality checks at run time
function sorttaxaCFperm!(pcf::Vector{Int8}, ptax::Vector{Int8})
    if ptax[1]==1
        if     ptax[2]==2
            pcf[1]=1
            if  ptax[3]==3 # ptax = 1,2,3,4
                pcf[2]=2; pcf[3]=3
            else           # ptax = 1,2,4,3
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==3
            pcf[1]=2
            if  ptax[3]==2 # ptax = 1,3,2,4
                pcf[2]=1; pcf[3]=3
            else           # ptax = 1,3,4,2
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==4
            pcf[1]=3
            if  ptax[3]==2 # ptax = 1,4,2,3
                pcf[2]=1; pcf[3]=2
            else           # ptax = 1,4,3,2
                pcf[2]=2; pcf[3]=1
            end
        end
    elseif ptax[1]==2
        if     ptax[2]==1
            pcf[1]=1
            if  ptax[3]==4 # ptax = 2,1,4,3
                pcf[2]=2; pcf[3]=3
            else           # ptax = 2,1,3,4
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==4
            pcf[1]=2
            if  ptax[3]==1 # ptax = 2,4,1,3
                pcf[2]=1; pcf[3]=3
            else           # ptax = 2,4,3,1
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==3
            pcf[1]=3
            if  ptax[3]==1 # ptax = 2,3,1,4
                pcf[2]=1; pcf[3]=2
            else           # ptax = 2,3,4,1
                pcf[2]=2; pcf[3]=1
            end
        end
    elseif ptax[1]==3
        if     ptax[2]==4
            pcf[1]=1
            if  ptax[3]==1 # ptax = 3,4,1,2
                pcf[2]=2; pcf[3]=3
            else           # ptax = 3,4,2,1
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==1
            pcf[1]=2
            if  ptax[3]==4 # ptax = 3,1,4,2
                pcf[2]=1; pcf[3]=3
            else           # ptax = 3,1,2,4
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==2
            pcf[1]=3
            if  ptax[3]==4 # ptax = 3,2,4,1
                pcf[2]=1; pcf[3]=2
            else           # ptax = 3,2,1,4
                pcf[2]=2; pcf[3]=1
            end
        end
    else # ptax[1]==4
        if     ptax[2]==3
            pcf[1]=1
            if  ptax[3]==2 # ptax = 4,3,2,1
                pcf[2]=2; pcf[3]=3
            else           # ptax = 4,3,1,2
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==2
            pcf[1]=2
            if  ptax[3]==3 # ptax = 4,2,3,1
                pcf[2]=1; pcf[3]=3
            else           # ptax = 4,2,1,3
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==1
            pcf[1]=3
            if  ptax[3]==3 # ptax = 4,1,3,2
                pcf[2]=1; pcf[3]=2
            else           # ptax = 4,1,2,3
                pcf[2]=2; pcf[3]=1
            end
        end
    end
end

"""
    setlengths!(edges::Vector{Edge}, lengths::Vector{Float64})

Assign new lengths to a vector of `edges`.
"""
@inline function setlengths!(edges::Vector{Edge}, lengths::Vector{Float64})
    for (e,l) in zip(edges, lengths)
        e.length = l
    end
end

"""
    getlengths(edges::Vector{Edge})

Vector of edge lengths for a vector of `edges`.
"""
getlengths(edges::Vector{Edge}) = [e.length for e in edges]

"""
    hashybridladder(net::HybridNetwork)

Return true if `net` contains a hybrid ladder: where a hybrid node's
child is itself a hybrid node.
This makes the network not treechild, assuming it is fully resolved.
(One of the nodes does not have any tree-node child).
"""
function hashybridladder(net::HybridNetwork)
    for h in net.hybrid
        if any(n.hybrid for n in getparents(h))
            return true
        end
    end
    return false
end

"""
    shrinkedge!(net::HybridNetwork, edge::Edge)

Delete `edge` from net, provided that it is a non-external tree edge.
Specifically: delete its child node (as determined by `isChild1`) and connect
all edges formerly incident to this child node to the parent node of `edge`,
thus creating a new polytomy, unless the child was of degree 2.

Warning: it's best for `isChild1` to be in sync with the root for this. If not,
the shrinking may fail (if `edge` is a tree edge but its "child" is a hybrid)
or the root may change arbitrarily (if the child of `edge` is the root).

Output: true if the remaining node (parent of `edge`) becomes a hybrid node with
more than 1 child after the shrinking; false otherwise (e.g. no polytomy was
created, or the new polytomy is below a tree node)
"""
function shrinkedge!(net::HybridNetwork, edge2shrink::Edge)
    edge2shrink.hybrid && error("cannot shrink hybrid edge number $(edge2shrink.number)")
    cn = getchild(edge2shrink)
    cn.hybrid && error("cannot shrink tree edge number $(edge2shrink.number): its child node is a hybrid. run directEdges! ?")
    pn = getparent(edge2shrink)
    isexternal(edge2shrink) &&  # (isleaf(cn) || isleaf(pn)) &&
      error("won't shrink edge number $(edge2shrink.number): it is incident to a leaf")
    removeEdge!(pn,edge2shrink)
    empty!(edge2shrink.node) # should help gc
    for ee in cn.edge
        ee !== edge2shrink || continue
        cn_index = findfirst(x -> x === cn, ee.node)
        ee.node[cn_index] = pn # ee.isChild1 remains synchronized
        push!(pn.edge, ee)
    end
    pn.hasHybEdge = any(e -> e.hybrid, pn.edge)
    empty!(cn.edge) # should help to garbage-collect cn
    deleteEdge!(net, edge2shrink; part=false)
    deleteNode!(net, cn)
    badpolytomy = false
    if pn.hybrid # count the number of pn's children, without relying on isChild1 of tree edges
        nc = sum((!e.hybrid || getchild(e) !== pn) for e in pn.edge)
        badpolytomy = (nc > 1)
    end
    return badpolytomy
end

@doc raw"""
    shrink2cycles!(net::HybridNetwork, unroot=false::Bool)

If `net` contains a 2-cycle, collapse the cycle into one edge of length
tA + γt1+(1-γ)t2 + tB (see below), and return true.
Return false otherwise.
A 2-cycle is a set of 2 parallel hybrid edges, from the same parent node to the
same hybrid child node.

           A                A
           | tA             |
         parent             |
           | \              |
    t2,1-γ |  | t1,γ        | tA + γ*t1 + (1-γ)*t2 + tB
           | /              |
         hybrid             |
           | tB             |
           B                B

If any of the lengths or gammas associated with a 2-cycle are missing,
the combined length is missing. If γ is missing, branch lengths
are calculated using γ=0.5.

If `unroot` is false and the root is up for deletion, it will be kept only if it
is has degree 2 or more. If `unroot` is true and the root is up for deletion, it
will be kept only if it has degree 3 or more. A root node with degree 1 will be
deleted in both cases.
"""
function shrink2cycles!(net::HybridNetwork, unroot=false::Bool)
    foundcycle = false
    nh = length(net.hybrid)
    ih = nh # hybrids deleted from the end
    while ih > 0
        h = net.hybrid[ih]
        minor = getparentedgeminor(h)
        major = getparentedge(h)
        pmin = getparent(minor) # minor parent node
        pmaj = getparent(major) # major parent node
        if pmin !== pmaj # no 2-cycle
            ih -= 1
            continue
        end
        # 2-cycle
        foundcycle = true
        shrink2cycleat!(net, minor, major, unroot)
        nh = length(net.hybrid)
        ih = nh
        # we re-do if a cycle was removed: a new cycle might have appeared
    end
    return foundcycle
end

"""
    shrink2cycleat!(net::HybridNetwork, minor::Edge, major::Edge, unroot::Bool)

Remove `minor` edge then update the branch length of the remaining `major` edge.
Called by [`shrink2cycles!`](@ref)

Assumption: `minor` and `major` do form a 2-cycle. That is, they start and end
at the same node.
"""
function shrink2cycleat!(net::HybridNetwork, minor::Edge, major::Edge,
                         unroot::Bool)
    g = minor.gamma
    if g == -1.0 g=.5; end
    major.length = addBL(multiplygammas(    g, minor.length),
                         multiplygammas(1.0-g, major.length))
    deletehybridedge!(net, minor, false,unroot,false,false,false) # nofuse,unroot,multgammas,simplify
    return nothing
end

"""
    shrink3cycles!(net::HybridNetwork, unroot=false::Bool)

Remove all 2- and 3-cycles from a network.

Return true if `net` contains a 2-cycle or a 3-cycle; false otherwise.
A 3-cycle (2-cycle) is a set of 3 (2) nodes that are all connected.
One of them must be a hybrid node, since `net` is a DAG.

If `unroot` is false and the root is up for deletion, it will be kept only if it
is has degree 2 or more. If `unroot` is true and the root is up for deletion, it
will be kept only if it has degree 3 or more. A root node with degree 1 will be
deleted in both cases.

See [`shrink3cycleat!`](@ref) for details on branch lengths and
inheritance values.
"""
function shrink3cycles!(net::HybridNetwork, unroot=false::Bool)
    foundcycle = false
    nh = length(net.hybrid)
    ih = nh # hybrids deleted from the end
    while ih > 0
        h = net.hybrid[ih]
        minor = getparentedgeminor(h)
        major = getparentedge(h)
        pmin = getparent(minor) # minor parent node
        pmaj = getparent(major) # major parent node
        if pmin === pmaj # 2-cycle
            foundcycle = true
            shrink2cycleat!(net, minor, major, unroot)
            nh = length(net.hybrid)
            ih = nh + 1 # start over if a cycle was removed by setting ih = nh + 1.
                        # Shrinking could have created a new cycle.
        else # 3-cycle
            result = shrink3cycleat!(net, h, minor, major, pmin, pmaj, unroot)
            if result
                foundcycle = true
                nh = length(net.hybrid)
                ih = nh + 1 # start over as above
            end
        end
        ih -= 1
    end
    return foundcycle
end

@doc raw"""
    shrink3cycleat!(net::HybridNetwork, hybrid::Node, edge1::Edge, edge2::Edge,
                    node1::Node, node2::Node, unroot::Bool)

Replace a 3-cycle at a given `hybrid` node by a single node, if any.
Assumption: `edge1` (`node1`) and `edge2` (`node2`) are the parent edges (nodes)
of `hybrid`. Return true if a 3-cycle is found and removed, false otherwise.
There is a 3-cycle if nodes 1 & 2 are connected, by an edge called `e3` below.

There are two cases, with differing effects on the γ inheritance
values and branch lengths.

**Hybrid case**: the 3-cycle is shrunk to a hybrid node, which occurs if
either node 1 or 2 is a hybrid node (that is, e3 is hybrid). If e3 goes
from node 1 to node 2, the 3-cycle (left) is shrunk as on the right:

    \eA      /eB           \eA  /eB
     1--e3->2       γ1+γ2γ3 \  / γ2(1-γ3)
      \    /               hybrid
     γ1\  /γ2
      hybrid

with new branch lengths:
new tA = tA + (γ1.t1 + γ2γ3.(t2+t3))/(γ1+γ2γ3),
new tB = tB + t2,
provided that γ1, γ2=1-γ1, and γ3 are not missing. If one of them is missing
then γ1 and γ2 remain as is, and e3 is deleted naively,
such that new tA = tA + t1 and new tB = tB + t2.
If γ's are not missing but one of t1,t2,t3 is missing, then the γ's are
updated to γ1+γ2γ3 and γ2(1-γ3), but t's are update naively.

**Tree case**: the 3-cycle is shrunk to a tree node, which occurs if node 1 & 2
are both tree nodes (that is, e3 is a tree edge). If eC is the child edge of
`hybrid`, the 3-cycle (left) is shrunk as on the right:

    \eA                  \eA
     1--e3--2--eB--       \
      \    /               n--eB--
     γ1\  /γ2              |
      hybrid               |eC
        |
        |eC

with new branch lengths:
new tA = tA + γ2.t3,
new tB = tB + γ1.t3,
new tC = tC + γ1.t1 + γ2.t2,
provided that γ1, γ2=1-γ1, t1, t2 and t3 are not missing.
If one is missing, then e1 is deleted naively such that
tB is unchanged, new tC = tC + t2 and new tA = tA + t3.
"""
function shrink3cycleat!(net::HybridNetwork, hybrid::Node, edge1::Edge,
                         edge2::Edge, node1::Node, node2::Node, unroot::Bool)
    # check for presence of 3 cycle
    edge3 = nothing
    for e in node1.edge # find edge connecting node1 and node2
        e !== edge1 || continue
        n = getOtherNode(e, node1)
        if n === node2
            edge3 = e
            break
        end
    end
    !isnothing(edge3) || return false # no 3-cycle at node h
    # identify case type
    if edge3.hybrid # one of the parent nodes is a hybrid
        # to shrink this, delete edge connecting these two nodes (edge3 here)
        if getchild(edge3) === node1
            node1, node2 = node2, node1
            edge1, edge2 = edge2, edge1
        end # now: node1 --edge3--> node2
        edgeA = nothing
        for e in node1.edge
            if e!== edge1 && e !== edge2
                edgeA = e
                break
            end
        end
        g1 = edge1.gamma
        g2g3 = multiplygammas(edge2.gamma, edge3.gamma)
        g1tilde = addBL(g1, g2g3)
        if g1tilde != -1.0 # else one of the γ is missing: do nothing with γs and ts
            edge1.gamma = g1tilde
            edge2.gamma = 1.0-g1tilde
            edge1.isMajor = g1tilde >= 0.5
            edge2.isMajor = !edge1.isMajor
            if edge1.length != -1.0 && edge2.length != -1.0 && edge3.length != -1.0
                edge1.length = (edge1.length *g1 + (edge3.length + edge2.length)*g2g3)/g1tilde
            end
        end
        deletehybridedge!(net, edge3, false,unroot,false,false) # nofuse,unroot,multgammas,simplify
    else # parent nodes 1 and 2 are both tree nodes
        edgeB = nothing
        for e in node2.edge
            if e !== edge1 && e !==edge3
                edgeB = e
                break
            end
        end
        g1t1 = multiplygammas(edge1.gamma, edge1.length)
        t3 = edge3.length
        if g1t1 != -1.0 && edge2.length != -1.0 && t3 != -1.0 # else do nothing: keep tA, tB, tC as is
            edgeB.length = addBL(edgeB.length, edge1.gamma * t3)
            edge3.length = t3 * edge2.gamma
            edge2.length = g1t1 + edge2.gamma * edge2.length
        end
        deletehybridedge!(net, edge1, false,unroot,false,false) # nofuse,unroot,multgammas,simplify
    end
    return true
end

"""
    adjacentedges(centeredge::Edge)

Vector of all edges that share a node with `centeredge`.
"""
function adjacentedges(centeredge::Edge)
    n = centeredge.node
    length(n) == 2 || error("center edge is connected to $(length(n)) nodes")
    @inbounds edges = copy(n[1].edge) # shallow copy, to avoid modifying the first node
    @inbounds for ei in n[2].edge
        ei === centeredge && continue # don't add the center edge again
        getOtherNode(ei, n[2]) === n[1] && continue # any parallel edge would have been in `edges` already
        push!(edges, ei)
    end
    return edges
end

#------------------------------------
function citation()
    bibfile = joinpath(@__DIR__, "..", "CITATION.bib")
    out = readlines(bibfile)
    println("Bibliography in bibtex format also in CITATION.bib")
    println(join(out,'\n'))
end
