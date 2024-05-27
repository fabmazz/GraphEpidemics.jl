using OffsetArrays
using Statistics

beta_R0(R0::Real,g::AbstractGraph,gamma::Real) = R0*gamma/mean(degree(g))
beta_R0(R0::Real,mean_deg::Real,gamma::Real) = R0*gamma/mean_deg

struct SIRSimData{I<:Integer, F<:AbstractFloat}
    rec_delays::Vector{I}
    infect_time::Vector{F}
    infect_node::Vector{F}
end


function explode(d::SIRSimData)
    (d.rec_delays, d.infect_time, d.infect_node)
end
function explode(d::SIRSimData, t::DataType)
    (convert.(t,d.rec_delays), convert.(t,d.infect_time), convert.(t,d.infect_node))
end

function SIRSimData(li::Union{Vector, Tuple})
    SIRSimData(li[1], li[2], li[3])
end

SIRSimData(rec_delays, N::Integer, typeDates::DataType) = SIRSimData(rec_delays, 
            fill(convert(typeDates,NaN),N), 
            fill(convert(typeDates,NaN),N))

function check_infection(rng::AbstractRNG, p::AbstractFloat, i::Integer, j::Integer, infect_prob_I::Symbol)
    rand(rng) < p
end

function check_infection(rng::AbstractRNG, p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)
    if infect_prob_IS == :I
        return rand(rng) < p[i_I]
    elseif infect_prob_IS == :S
        return rand(rng) < p[j_S]
    else
        return rand(rng) < sqrt(p[i_I]*p[j_S])
    end
end

function get_p_infection(p::AbstractFloat, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)
    p
end

function get_p_infection(p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)
    if infect_prob_IS == :I
        return p[i_I]
    elseif infect_prob_IS == :S
        return p[j_S]
    else
        return sqrt(p[i_I]*p[j_S])
    end
end

function get_p_infection(p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::InfectI)
    return p[i_I]
end
function get_p_infection(p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::InfectS)
    return p[j_S]
end
function get_p_infection(p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::InfectMeanIS)
    return sqrt(p[j_S]*p[i_I])
end
function get_p_infection(p::AbstractFloat, i_I::Integer, j_S::Integer, infect_prob_IS::InfectionDirection)
    return p
end




function sim_sir_fast(g::AbstractGraph, model::SIRModel, T::Integer, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{I}; beta_IorS::InfectionDirection = InfectI()) where I<: Integer
    N = nv(g)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    delays = simdata.rec_delays

    #=if !(beta_IorS in (:I,:S,:mean))
        throw(ArgumentError("argument 'beta_IorS' must be one of :I, :S, or :mean"))
    end
    =#

    states::Vector{Int8} = fill(1, N)
    for i in patient_zeros
        states[i] = 2
        infect_t[i] = -1 ## time of infection
        infect_i[i] = -10 ## node of infection
    end
    
    
    t=0
    useT = (T >= 0)
    if useT
        counts = Vector{NTuple{3,Int}}(undef, T+1)
    else
        counts = Vector{NTuple{3,Int}}(undef, 0)
    end
    nI = length(patient_zeros)
    mcont = true
    while (mcont)
        ## set recovered state
        #println("Have $(sum(states.==1)) infected, $(sum(states.==0)) sus PRE")
        mask_rec = @. t >= ( infect_t + 1 + delays)
        states[mask_rec] .= 3
        nS =  (sum(states.==1))
        nI = sum(states.==2)
        nR= sum(mask_rec)
        if useT
            counts[t+1] = (nS, nI, nR)
            ## set here the flag
            mcont = (t+1 <= T)
        else 
            push!(counts, (nS, nI,nR))
            mcont = nI > 0
        end
        #counts[t+1,1] = nS
        #counts[t+1,2] = nI
        #counts[t+1,3] = nR
        #println("Have $nS S $(sum(states.==1)) I $nR R")
        c=0
        for i in findall(states.==2)
            ## TODO: parallel runs need to lock the graph object and then unlock it
            for j in neighbors(g,i)
                if (states[j] == 1) & isnan(infect_i[j]) ## double check to be sure
                    ## Here, multiple dispatch will help distinguish 
                    ## when I have a vector of probabilities or just a single one for everyone
                    if rand(rng) < get_p_infection(model.beta, i, j, beta_IorS)
                        ## infected
                        infect_t[j] = t
                        infect_i[j] = i
                        states[j] = 2
                        c+=1
                    end
                end
            end #for neighbors
        end #for infected
        #println("$c new infected at time $t")
        
        t+=1
    end #for time loop
    states, counts
end

function run_sir_fast(g::AbstractGraph, model::SIRModel, T::Integer, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; beta_IorS::InfectionDirection = InfectI(), dtype::DataType=Float64)
    ## draw recovery delays
    N= nv(g)
    nodes = collect(Int32,1:N)
    delays = draw_delays(model, model.gamma, rng, nodes)

    data = SIRSimData(delays,N, dtype)
    endstate, cc = sim_sir_fast(g, model, T, data, rng, patient_zeros, beta_IorS=beta_IorS)

    data, cc
end

function calc_n_comparts(infect_times,delays,T::Integer)
    ti= infect_times .+ 1
    isinfect= collect(0:T) .>= (ti)'
    isrec = collect(0:T) .>= (ti.+delays)'
    nR = sum(isrec,dims=2)
    nI = sum(isinfect,dims=2).-nR
    nS = sum(.~isinfect,dims=2)
    hcat(nS,nI,nR)
end

function calc_n_states(infects, delays, T::Integer)
    #counts = zeros(T+1) ## from 0 to T
    counts = OffsetArray(zeros(Int,T+1,3), 0:T,1:3)

    N = length(delays)
    for i=1:N
        if isnan(infects[i])
            ## never infected
            counts[:,1].+=1
            continue
        end
        ti = convert(Int,infects[i])+1
        if ti < 0
            ## never infected
            counts[:,1] .+= 1
            continue
        end
        trec = min(ti + delays[i], T)
        counts[0:(ti-1),1] .+=1 ## times S
        counts[ti:(trec-1),2] .+= 1
        counts[trec:end,3].+=1
    end
    counts
end

function calc_nstates_all(infects,recovs,T::Integer)
    dats = [calc_n_comparts(ii,dd, T) for (ii,dd) in zip(infects,recovs)]

    counts=cat(dats..., dims=3)
end

function calc_Rarr(inf,N, T, tmax, win_size=7)
    newinf = zeros(Int,T+2)
    succinf = zeros(Int,T+2)
    
    for i=1:N
        if isnan(inf[i,1])
            continue
        end
        ti=convert(Int,inf[i,1])
        newinf[ti+2]+=1
        if !isnan(inf[i,2]) && inf[i,2]>=0
            ji = convert(Int,inf[i,2])
            
            tj = convert(Int, inf[ji,1])
            if tj >= 0
                succinf[tj+2]+=1
            end
        end
    end
    Rarr = zeros(tmax)
    ws = win_size
    for t=ws:tmax
        Rarr[t]=sum(succinf[max(t-ws,1):t])/sum(newinf[max(t-ws,1):t]) 
    end
    newinf, succinf, Rarr
end