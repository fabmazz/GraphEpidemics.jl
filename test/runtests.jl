using DataFrames
using Random
using Distributions
using Graphs
using GraphEpidemics
using Test

import GraphEpidemics: model_states, spreading_states, trans_independent, draw_delays, first_active_states, draw_delays_markov

struct SEIRModel{F<:AbstractFloat} <: AbstractEpiModel
    eta::F
    beta::F ##vector is the 
    gamma::F
    #stateType::DataType
end

model_states(x::SEIRModel) = (:S,:E,:I,:R)
spreading_states(x::SEIRModel) = Dict(:I=>[(:S,:E, x.beta)])
trans_independent(x::SEIRModel) = [(:E,:I,x.eta),(:I,:R, x.gamma)]

first_active_states(m::SEIRModel) = (:E,)
draw_delays(m::SEIRModel, p, rng::AbstractRNG, node) = draw_delays_markov(p, rng, node )

function test_run_SEIR()
    model=SEIRModel(1/5, 0.06, 0.11)

    rng = Xoshiro(40)
    g=erdos_renyi(2000,0.004, rng=rng)
    T=140
    is_connected(g)
    
    data = init_model_discrete(model, g,rng,:S)
    start=rand(rng,vertices(g),2)

    set_state_nodes(model,data,rng,start,:E,true)
    counts = (run_complex_contagion(model,g,T,rng, data));

    counts[end] == (; :t=>140, :I=>3, :R=>1891, :S=>106,:E=>0)
end
@testset "SEIR" begin
    
    @test test_run_SEIR()


end

@testset "recurrent" begin
    model=SIRSModel(0.1,0.1,0.05)

    rng = Xoshiro(40)
    g=barabasi_albert(1000,14, rng=rng)

    #println(is_connected(g))

    #@benchmark  begin
    data = init_model_discrete(model, g,rng,:S)
    start=rand(rng,vertices(g),2)

    set_state_nodes(model,data,rng,start,:I,true)
    counts = (run_complex_contagion(model,g,200,rng, data));

    count_df = DataFrame(counts)
    l = size(count_df,1)
    r=count_df[(end-50):end,:]
    select!(r,Not(:t))

    @test round.(Int,mean.(eachcol(r))) == [316, 618,66]
end