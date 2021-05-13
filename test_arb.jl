using JuMP
using Gurobi
using MAT
using Printf
# using Plots
using JLD

# High-level Settings
Zone = "NYC" # price zone name

# read price
fileln = matopen(string("./RTP_data/RTP_",Zone,"_2010_2019.mat"))
RTP = read(fileln, "RTP")
close(fileln)

#

E = 1;  # storage energy capacity in MWh
# simulation setting
T = 288; # time step per day
M = 1/12; # duaration per step in hour
N_sim = 3650; # number of days
L = RTP[:,1];

Lp = L .> 0

# BES setting
P = .25; # power rating MW
# E = 2; # energy rating MWh
eta = .9; # single-trip efficiency
e0 = .5 * E;
ef = e0;
MC = 0; # marginal discharge cost

# initialize optimization model
model = Model(Gurobi.Optimizer)
set_silent(model) # no outputs
# total discharge power
@variable(model, d[1:T], lower_bound = 0)
# total charge power
@variable(model, c[1:T], lower_bound = 0)
# total energy level
@variable(model, e[1:T], lower_bound = 0)
@variable(model, D) # incremental cycle degradation
@variable(model, C) # value of battery capacity at the end of operation
@variable(model, R) # market revenue

# arbitrage revenue
@constraint(model, ArbRev, R == M*sum(L.*(d-c)) )
# degradation increments
# piece-wise linear degradation opportunity value
@constraint(model, DegCost, C == M*MC*sum(d)  )
# initial SoC evolution
@constraint(model, SoCInit, e[1] - e0 == M*(c[1] - d[1]) )
# rest SoC evolution
@constraint(model, SoCCont[t = 2:T], e[t] - e[t-1] == M*(c[t]*eta - d[t]/eta) )
# final energy level
@constraint(model, Enelast, e[T] >= ef )
# max discharge power
@constraint(model, DisMax[t=1:T], d[t] <= P*Lp[t] )
# max charge power
@constraint(model, ChrMax[t=1:T], c[t] <= P )
# max energy level
@constraint(model, SoCMax[t=1:T], e[t] <= E )
# warranty constraint
# @constraint(model, ETHMax, M*sum(d)/eta <= E-E*D0 )

# maximize revenue plus degradation value
@objective(model, Max, R-C)



# initialize
R_s = zeros(1, N_sim)
P_s = zeros(1, N_sim)
C_s = zeros(1, N_sim)
D_s = zeros(1, N_sim)

@printf("Optimization starts...\n")
for n = (N_sim-364):(N_sim)


# update prices
local L = RTP[:,n]
local Lp = L .> 0

# update prices in constraints
for t = 1:T
    set_normalized_coefficient(ArbRev, d[t], -M*L[t] )
    set_normalized_coefficient(ArbRev, c[t], M*L[t] )
    set_normalized_rhs(DisMax[t], P*Lp[t])
end


optimize!(model)

global R_s[n] = value(R)# objective_value(model);

global C_s[n] = value(C)
global P_s[n] = value(R)-value(C)
global D_s[n] = sum(value.(d))


termination_status(model)
@printf("Finished Day %d, Cum Rev %d, Cum Profit %d, Cum Cost %d, OptStatus: %s \n", n, sum(R_s), sum(P_s), sum(C_s), termination_status(model))


end


# FN = @sprintf("H%d_L%d_P%d_Ben", E, EoL*100, P_recycle/1e4)
#
#
# file = matopen(string("./results_mat/",Zone,"2010_2019_",FN,".mat"), "w")
# write(file, "D_s", collect(D_s))
# write(file, "E", E)
# write(file, "P", P)
# write(file, "eta", eta)
# write(file, "J", J)
# write(file, "Ycal", Ycal)
# write(file, "P_recycle", P_recycle)
# write(file, "Vs0", Vs0)
# write(file, "EoL", EoL)
# write(file, "EoW", EoW)
# close(file)

# save(string(pwd(), "\\results_jl\\",Zone,"2010_2019_",FN,".jld"), "Vs", Vs, "Ds", Ds, "P", P, "E", E)



# end
