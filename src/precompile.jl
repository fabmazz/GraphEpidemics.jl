using StatsBase: sample
N=100
d=10
G=random_regular_graph(N,d, seed=4)

DEG=d

rng = Xoshiro(2)
gam = 1/10
bet = beta_R0(1.8,DEG,gam)
a1 = rand(rng,N)*0.3 .+0.7
s1 = rand(rng,N)*0.5 .+0.5
model = SEIRHetModel(bet,a1,s1,gam*2,gam )

sd,cc = run_seir_fast(G,model,-1,rng, sample(rng,1:N,5,replace=false),delays_type=Int32,dtype=Float32)
cpl = stack(cc)



model = SIRModel(bet,gam)

run_sir_fast(G, model, -1, rng, sample(rng,1:N,5,replace=false))