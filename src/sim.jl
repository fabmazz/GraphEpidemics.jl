using OffsetArrays

beta_R0(R0::Real,g::AbstractGraph,gamma::Real) = R0*gamma/mean(degree(g))
beta_R0(R0::Real,mean_deg::Real,gamma::Real) = R0*gamma/mean_deg


function sim_sir(g::AbstractGraph,T::Integer, lambda::AbstractFloat, delays::Vector, rng::AbstractRNG, pz)
    N = nv(g)
    #track = fill(NaN, (N,2))
    infect_t = fill(NaN, N)
    infect_i = fill(NaN, N)

    states::Vector{Int8} = fill(0, N)
    for i in pz
        states[i] = 1
        infect_t[i] = -1 ## time of infection
        infect_i[i] = -10 ## node of infection
    end

    for t =0:T
        ## set recovered state
        #println("Have $(sum(states.==1)) infected, $(sum(states.==0)) sus PRE")
        #states[@. t >= ( track[:,1] +1)] .= 1
        states[@. t >= ( infect_t + 1 + delays) ] .= 2
        nS =  (sum(states.==0))
        nR= sum(states.==2)
        #println("Have $nS S $(sum(states.==1)) I $nR R")
        c=0
        for i in findall(states.==1)
            for j in neighbors(g,i)
                if (states[j] == 0) & isnan(infect_i[j]) ## double check to be sure
                    if rand(rng) < lambda
                        ## infected
                        infect_t[j] = t
                        infect_i[j] = i
                        states[j] = 1
                        c+=1
                    end
                end
            end #for neighbors
        end #for infected
        #println("$c new infected at time $t")

    end #for time loop
    infect_t, infect_i, states
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