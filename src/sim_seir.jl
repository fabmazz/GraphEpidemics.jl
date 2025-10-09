
struct SEIRSimData{I<:Real, F<:AbstractFloat}
    lat_delays::Vector{I}
    rec_delays::Vector{I}
    infect_time::Vector{F}
    infect_node::Vector{F}
end
SEIRSimData{I,F}(N::Integer) where {I,F}= SEIRSimData(zeros(I,N),zeros(I,N), fill(F(NaN),N),fill(F(NaN),N))
SEIRSimData{F}(lat_delays::Vector{I}, rec_delays::Vector{I}, N::Integer) where {I, F<:AF} = SEIRSimData(
    lat_delays, rec_delays, fill(F(NaN), N), fill(F(NaN,N))
) 

function reset!(data::SEIRSimData, model::AbstractSEIRModel)
    draw_delay_geom!(gamma(model), rng, data.rec_delays)
    draw_delay_geom!(eps(model), rng, data.lat_delays)
    fill!(data.infect_time, NaN)
    fill!(data.infect_node, NaN)
end
"""
`SEIRSimData(lat_delays, rec_delays, N::Integer, typeDates::DataType)`

Constructs a `SEIRSimData` object given latent delays, recovery delays, number of nodes, using the specified time data type.
DEPRECATED!
"""
SEIRSimData(lat_delays, rec_delays, N::Integer, typeDates::DataType) = SEIRSimData(lat_delays,rec_delays, 
            fill(convert(typeDates,NaN),N), 
            fill(convert(typeDates,NaN),N))

struct SEIRStatesCounter <:AbstractStatesCounter
end

function count_states(counter::SEIRStatesCounter,states::Vector)
    c=SVector{4,Int}(
        sum(states.==i) for i=1:4
    )
end

"""
`sim_seir_fast(g::AbstractGraph, model::AbstractSIRModel, T::Integer, simdata::SIRSimData, rng::AbstractRNG, patient_zeros::Vector{I}; beta_IorS=:I, counter=BaseSIRStatesCounter()) where I<:Integer`

Runs a discrete-time SEIR simulation over `T` time steps, with the data structures already set up.

### Optional arguments
- `counter=BaseSIRStatesCounter()`: The state-counting struct. When using a custom struct, make sure to implement `count_states(::BaseSIRStatesCounter, states::Vector)`
- `dynstateChanger`: object used for changing the state dynamically
"""
function sim_seir_fast(G::AbstractGraph, model::AbstractSEIRModel, T::Integer, simdata::SEIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{I}; 
    counter::AbstractStatesCounter=SEIRStatesCounter(), 
    dynstateChanger::Union{Nothing,AbstractStateChanger}=nothing) where I<: Integer
    N = nv(G)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    del_rec = simdata.rec_delays
    del_lat = simdata.lat_delays

    states::Vector{StI} = fill(1, N)
    for i in patient_zeros
        states[i] = 2
        infect_t[i] = -1 ## time of infection
        infect_i[i] = -10 ## node of infection
    end
    
    
    t=0
    useT = (T >= 0)
    if useT
        trace_states = Vector{SVector}(undef, T+1)
    else
        trace_states = Vector{SVector}(undef, 0)
    end
    nE = length(patient_zeros)
    Eidx = Set(patient_zeros)
    mcont = true
    nI = 0#length(patient_zeros)
    Iidx = Set()
    new_expo = Set()
    while (mcont)
        
        new_inf = Set()
        for i in Eidx
            if (t>=infect_t[i]+1+del_lat[i])
                states[i] = 3
                push!(new_inf, i)
            end
        end
        new_rec = Set()
        for i in Iidx
            if (t >= infect_t[i] + del_lat[i] + del_rec[i]+1)
                states[i] = 4
                push!(new_rec, i)
            end
        end
        #states[mask_rec] .= 3
        #nS =  (sum(states.==1))
        setdiff!(Iidx, new_rec)
        setdiff!(Eidx, new_inf)
        union!(Iidx, new_inf)
        union!(Eidx, new_expo)
        nI = sum(states.==3)
        nE = length(Eidx)
        @assert nI  == length(Iidx)
        #nR= N - nI - nS 
        if useT
            # use t+1 because of indexing
            trace_states[t+1] = count_states(counter,states) #[nS, nI, nR]
            ## set here the flag
            mcont = (t+1 <= T)
        else 
            push!(trace_states, count_states(counter,states))
            mcont = (nI > 0) || (nE > 0)
        end

        ### State changer
        if (!isnothing(dynstateChanger))
            change_states_dyn(dynstateChanger, model, states, G, t, Iidx)
        end

        c=0
        #Iidx = findall(states.==2) # this will be used also later
        empty!(new_expo)
        for i in Iidx
            ## TODO: parallel runs need to lock the graph object and then unlock it
            for j in neighbors(G,i)
                if (states[j] == 1) ## double check to be sure
                    ## Here, multiple dispatch will help distinguish 
                    ## when it's a vector of probabilities or just a single one for everyone
                    if rand(rng) < calc_prob_infection(model, i, j, :ignored)
                        ## infected
                        infect_t[j] = t
                        infect_i[j] = i
                        states[j] = 2
                        push!(new_expo, j)
                        c+=1
                    end
                end
            end #for neighbors
        end #for infected
        #println("$c new infected at time $t")
        
        t+=1
    end #for time loop
    states, trace_states
end
"""
`run_seir_fast(g::AbstractGraph, model::AbstractSIRModel, T::Integer, rng::AbstractRNG, patient_zeros::Vector{<:Integer}; dtype=Float64,  delays_types::DataType=Int64, counter=BaseSIRStatesCounter())`

Sets up and runs a fast SEIR and SEIR-like simulation, returning simulation data and state counts.

### Optional arguments
- `dtype=Float64`: Data type for output arrays.
- `delays_types=Int64`: Data type for the delays
- `counter=BaseSIRStatesCounter()`: State counting method.
- `dynstateChanger`: object used for changing the state dynamically

"""
function run_seir_fast(g::AbstractGraph, model::AbstractSEIRModel, T::Integer, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; dtype::DataType=Float64, delays_type::DataType=Int64,
    counter::AbstractStatesCounter=SEIRStatesCounter(), dynstateChanger::Union{Nothing, AbstractStateChanger}=nothing)
    ## draw delays
    N= nv(g)
    # nodes = collect(Int32,1:N)
    rec_delays = Vector{delays_type}(undef,N)
    draw_delay_geom!(gamma(model), rng, rec_delays)
    latent_delays = Vector{delays_type}(undef,N)
    draw_delay_geom!(eps(model), rng, latent_delays)

    #draw_delays(model, gamma(model), rng, nodes)
    #latent_delays = draw_delays(model, eps(model), rng,nodes)

    data = SEIRSimData(latent_delays,rec_delays,N, dtype)
    endstate, cc = sim_seir_fast(g, model, T, data, rng, patient_zeros,counter=counter, dynstateChanger=dynstateChanger)

    data, cc
end