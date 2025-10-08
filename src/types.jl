#= 
 Copyright (c) 2025 Fabio Mazza
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at https://mozilla.org/MPL/2.0/.
=#

"""
`AbstractEpiModel`

Abstract supertype for all epidemiological models. To be used in `run_complex_contagion`, subtypes must implement
- `model_states(::model)` which returns a tuple of symbols that describe the model states
- `spreading_states(::model)` returning a dict of transitions describing the infector state, the state of soon to be infected (from and to) and the rate, e.g. `Dict(:I=>[(:S,:I, x.beta)])`
- `trans_independent(::model)` which is a list of all the independent transitions with their rate (e.g., `[(:I,:R, x.gamma)]` for the SIR model
- `first_active_states(::model)` which describes the seed states (e.g. `(:I,)` for the SIR)
"""
abstract type AbstractEpiModel end

"""
`AbstractSIRModel`

Abstract supertype for all SIR epidemiological models. Subtypes must implement a method for computing infection probabilities.
This model is used for `run_sir_fast` and `run_sir_gillespie`.
"""
abstract type AbstractMarkovModel <: AbstractEpiModel end

abstract type AbstractSIRModel <: AbstractMarkovModel end
abstract type AbstractSEIRModel <: AbstractMarkovModel end


abstract type AbstractStatesCounter end
abstract type AbstractStateChanger end

struct BaseSIRStatesCounter <:AbstractStatesCounter
end

function count_states(counter::BaseSIRStatesCounter,states::Vector)
    c=SVector{3,Int}(
        sum(states.==i) for i=1:3
    )
end

function change_states_dyn(dynstateChanger::AbstractStateChanger, model::AbstractSIRModel, states, g, t, Iidx)
    ### EMPTY
end

struct NotImplementedError <: Exception
    message::String
end

Base.showerror(io::IO, e::NotImplementedError) = print(io, e.message)

const AF=AbstractFloat
