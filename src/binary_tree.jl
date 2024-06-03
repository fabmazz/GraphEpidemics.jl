height_tree(n_leaves::Real) = ceil(Int, log2(n_leaves)) +1 

function calc_n_elemns(n::Integer)
   nl = ceil(Int,log2(n))+1
   2^nl-1
end
    
function find_h_rec(idc::Integer, lvlmax::Integer)
    #println("idc $idc lvl $lvlmax")
    if idc == 3
        return 2
    elseif idc == 2
        return 1
    elseif idc == 1
        return 1
    elseif idc <= 0
        return -10
    end

    if idc > 2^lvlmax-1
        throw(DomainError((idc,lvlmax),"The index must be lower than the power specified by the maximum height"))
    elseif idc == 2^lvlmax -1
        return lvlmax
    elseif idc > 2^(lvlmax-1)-1
        return find_h_rec(idc- 2^(lvlmax-1)+1, lvlmax-1)
    else
        return find_h_rec(idc, lvlmax-1)
    end
end

function idcs_leaves(lvl::I) where I <: Integer
    if lvl <=1 
        return I[]
    end
    idcs = I[1,2]
    l = 2
    while l < lvl
        idcs = vcat(idcs, idcs.+ (2^l -1))
        l+=1
    end
    idcs
end

function up(i::Integer, l::Integer, isleft::Bool)
    if isleft
        return i+2^l
    else
        return i+1
    end
end

function down(i::Integer, l::Integer, left::Bool)
    if left
        return i- 2^(l-1)
    else
        return i-1
    end
end

function calc_isleft_vector(height::Integer)
    if height ==1
        return [true]
    end
    v = [true,false,true]
    h = 2
    while (h<height)
        v = vcat(v,v)
        v[end] = false
        push!(v,true)
        h+=1
    end
    v
end

#=
struct Event1
    idx::Int
    from::Int8
    to::Int8
    delay::Float64
end

Event=Event1
rate(e::Event1) = e.delay
hash(e::Event1) = hash((e.idx,e.from,e.to))
isequal(e1::Event, e2::Event) = (e1.idx==e2.idx) && (e1.from==e2.from) && (e1.to==e2.to)
=#

@Base.kwdef mutable struct BinaryTree{F<:AbstractFloat, T}
    weights::Vector{F}
    idxfree::Set{Int}
    noccupied::Int
    height::Int
    numleaves::Int
    leftchild::Vector{Bool}
    eventsByPos::Dict{Int, T}
    posOfEvent::Dict{T,Int}
    function BinaryTree{F,T}(nleaves::Int) where {F,T}
        h = height_tree(nleaves)
        nels = calc_n_elemns(nleaves)
        ileaves = idcs_leaves(h)
        new{F,T}(zeros(nels), Set(ileaves), 0, h, 2^h, calc_isleft_vector(h), Dict{Int, T}(), Dict{T, Int}())
    end
end

function BinaryTree(v::Vector{T}) where T
    nleaves = max(length(v),2)
    h = height_tree(nleaves)
    nels = calc_n_elemns(nleaves)
    ileaves = idcs_leaves(h)
    tree=BinaryTree(zeros(nels), Set(ileaves), 0, h, 2^h, calc_isleft_vector(h), Dict{Int, T}(), Dict{T, Int}())
    for e in v
        add_event!(tree,e)
    end
    tree
end

function expand!(tree::BinaryTree)
    h = tree.height
    # new empty indices
    new_empties = (idcs_leaves(h).+(2^h -1 ))
    union!(tree.idxfree, new_empties)
    
    current_sum = tree.weights[end]
    tree.weights= vcat(tree.weights, zeros(length(tree.weights)))
     #append!(tree.weights, zeros(length(tree.weights)))
    push!(tree.weights, current_sum)
    
    new_left = copy(tree.leftchild)
    new_left[end] = false
    append!(tree.leftchild, new_left)
    push!(tree.leftchild, true)

    tree.height += 1
    tree.numleaves = 2^(tree.height)
    tree
end

function update_weights!(tree::BinaryTree, imodified::Integer, lvlmod::Integer)
    if imodified == length(tree.weights)
        return
    end
    ni = up(imodified,lvlmod, tree.leftchild[imodified])
    li = down(ni, lvlmod+1, !tree.leftchild[imodified])
    tree.weights[ni] = tree.weights[li]+tree.weights[imodified]
    update_weights!(tree, ni, lvlmod+1)
end

function add_event!(tree::BinaryTree, ev, sum_tree::Bool=true) 
    if ev in keys(tree.posOfEvent)
        throw(ArgumentError("The event $ev is already in the tree. Current keys: $(keys(tree.posOfEvent))"))
    end
    if length(tree.idxfree) == 0
        expand!(tree)
    end

    ni = pop!(tree.idxfree)
    
    tree.weights[ni] = rate(ev)

    tree.posOfEvent[ev] = ni
    tree.eventsByPos[ni] = ev

    tree.noccupied+=1

    if sum_tree
        update_weights!(tree, ni, 1)
    end
end

function find_leaf_idx_random_draw(t::BinaryTree,v::Real,debug=false)
    istart = length(t.weights) #calc_n_elemns(t.noccupied)
    h = find_h_rec(istart, t.height)
    itop = istart
    while (h > 1)
        idl = down(itop, h, true)
        idr = down(itop, h, false)
        if v <= 0
            throw(DomainError(v, "The value to compare must be positive"))
        end
        if debug
            println("h $h, v: $v. The left value is $(t.weights[idl]), right is $(t.weights[idr])")
        end
        ## these two choices are to be checked first
        if t.weights[idl] == 0
            itop = idr
        elseif  t.weights[idr] == 0
            itop = idl
        #println("i $itop, h $h")
        elseif v <= t.weights[idl]
            itop = idl
        else
            v = v-t.weights[idl]
            itop = down(itop, h, false)
        end
        
        h -= 1
    end
    itop
end
    
function remove_event_idx!(t::BinaryTree, idx)
    ev = t.eventsByPos[idx]
    delete!(t.eventsByPos, idx)
    delete!(t.posOfEvent, ev)
    push!(t.idxfree, idx)

    t.weights[idx] = 0.0
    update_weights!(t, idx,1)
    t.noccupied -=1

    ev
end
function remove_event!(t::BinaryTree, event)
    idx = t.posOfEvent[event]
    remove_event_idx!(t,idx)
    idx
end

function update_event!(t::BinaryTree, event, r::AbstractFloat=-1.0)
    if !(event in keys(t.posOfEvent))
        throw(ArgumentError("the event $event is not in the dict, add it instead"))
    end
    idx = t.posOfEvent[event]
    t.weights[idx] = r>0 ? r : rate(event)
    update_weights!(t, idx, 1)

    t.eventsByPos[idx] = event
end


function draw_event_random(rng::AbstractRNG, tree::BinaryTree)
    imax = calc_n_elemns(t.noccupied)
    r = rand(rng) * tree.weights[imax]

    i_w = find_leaf_idx_random_draw(tree, r)

    ev = t.eventsByPos[i_w]

    ev
end

total_rate(t::BinaryTree) = t.weights[end]
is_empty(t::BinaryTree) = (t.noccupied == 0)
has_event(t::BinaryTree, ev ) = haskey(t.posOfEvent, ev)

get_event_by_idx(t::BinaryTree, idx) = t.eventsByPos[idx]