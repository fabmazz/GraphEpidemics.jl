# GraphEpidemics

*Run complex epidemic models on graphs*

## Overview

This package aims to provide a general method for running different epidemic models on complex graphs, in a very efficient way. 

The models run in discrete time intervals, with each node of the graph representing an individual that can be in a particular (user defined) state.

The graphs are provided via the [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) Julia package

## Installation

This package is not yet included in the Julia general index.

To install it, from the Julia Pkg manager (`]`), run
```
pkg> add https://github.com/fabmazz/GraphEpidemics.jl.git
```
or, if you want to contribute to the development, first checkout the repository, and then navigate to the folder above and run 
```
pkg> dev ./GraphEpidemics.jl
```
which will keep the package source in the directory where you downloaded it.

## Usage

You can look at the [example notebook](https://github.com/fabmazz/GraphEpidemics.jl/blob/main/example/Models_example.ipynb) provided for an example of running a standard SIR Model and an implementation of the SEIR model.

## Caveats

At the moment, only non-recurring epidemic models are supported.
