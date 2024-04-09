using OffsetArrays
using Statistics

beta_R0(R0::Real,g::AbstractGraph,gamma::Real) = R0*gamma/mean(degree(g))
beta_R0(R0::Real,mean_deg::Real,gamma::Real) = R0*gamma/mean_deg

struct SIRSimData{I<:Integer, F<:AbstractFloat}
    rec_delays::Vector{I}
    infect_time::Vector{F}
    infect_node::Vector{F}
end

SIRSimData(rec_delays, N::Integer) = SIRSimData(rec_delays, fill(NaN,N), fill(NaN, N))

function check_infection(rng::AbstractRNG, p::AbstractFloat, i::Integer, j::Integer, infect_prob_I::Bool)
    rand(rng) < p
end

function check_infection(rng::AbstractRNG, p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_I::Bool)
    if infect_prob_I
        return rand(rng) < p[i_I]
    else
        return rand(rng) < p[j_S]
    end
end

function sim_sir_fast(g::AbstractGraph, model::SIRModel, T::Integer, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{I}; infect_prob_I::Bool = true) where I<: Integer
    N = nv(g)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    delays = simdata.rec_delays

    states::Vector{Int8} = fill(1, N)
    for i in patient_zeros
        states[i] = 2
        infect_t[i] = -1 ## time of infection
        infect_i[i] = -10 ## node of infection
    end

    counts = zeros(Int,(T+1,3))
    for t =0:T
        ## set recovered state
        #println("Have $(sum(states.==1)) infected, $(sum(states.==0)) sus PRE")
        mask_rec = @. t >= ( infect_t + 1 + delays)
        states[mask_rec] .= 3
        nS =  (sum(states.==1))
        nI = sum(states.==2)
        nR= sum(mask_rec)
        counts[t+1,1] = nS
        counts[t+1,2] = nI
        counts[t+1,3] = nR
        #println("Have $nS S $(sum(states.==1)) I $nR R")
        c=0
        for i in findall(states.==2)
            for j in neighbors(g,i)
                if (states[j] == 1) & isnan(infect_i[j]) ## double check to be sure
                    ## Here, multiple dispatch will help distinguish 
                    ## when I have a vector of probabilities or just a single one for everyone
                    if check_infection(rng, model.beta, i, j, infect_prob_I)
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

    end #for time loop
    states, counts
end

function run_sir_fast(g::AbstractGraph, model::SIRModel, T::Integer, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; prob_infect_I::Bool=true)
    ## draw recovery delays
    N= nv(g)
    nodes = collect(Int32,1:N)
    delays = draw_delays(model, model.gamma, rng, nodes)

    data = SIRSimData(delays,N)
    endstate, cc = sim_sir_fast(g, model, T, data, rng, patient_zeros, infect_prob_I=prob_infect_I)

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