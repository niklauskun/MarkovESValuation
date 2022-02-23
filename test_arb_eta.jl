using JuMP
using Gurobi
using MAT
using Printf
using DataFrames
using CSV
# using Plots
using JLD

# High-level Settings
Zone = "NYC" # price zone name

# read price
fileln = matopen(string("./RTP_data/RTP_",Zone,"_2010_2019.mat"))
RTP = read(fileln, "RTP")
close(fileln)
#

E1 = 0.2;  # storage energy capacity in MWh
E2 = 0.7;
E3 = 0.1;
# simulation setting
T = 288; # time step per day
M = 1/12; # duaration per step in hour
N_sim = 3650; # number of days
L = RTP[:,1];

Lp = L .> 0

# BES setting
P = .5; # power rating MW
# E = 2; # energy rating MWh
eta1 = .8; # single-trip efficiency
eta2 = .9;
eta3 = .7
e01 = .2;
e02 = .3;
e03 = 0;
ef = e01 + e02 + e03;
# ef = .8 * E;
MC = 10; # marginal discharge cost

# initialize optimization model
model = Model(Gurobi.Optimizer)
set_silent(model) # no outputs
# total discharge power
@variable(model, d[1:T], lower_bound = 0)
# total charge power
@variable(model, c[1:T], lower_bound = 0)
# E1 discharge power
@variable(model, d1[1:T], lower_bound = 0)
# E1 charge power
@variable(model, c1[1:T], lower_bound = 0)
# E1 energy level
@variable(model, e1[1:T], lower_bound = 0)
# E2 discharge power
@variable(model, d2[1:T], lower_bound = 0)
# E2 charge power
@variable(model, c2[1:T], lower_bound = 0)
# E2 energy level
@variable(model, e2[1:T], lower_bound = 0)
# E3 discharge power
@variable(model, d3[1:T], lower_bound = 0)
# E3 charge power
@variable(model, c3[1:T], lower_bound = 0)
# E3 energy level
@variable(model, e3[1:T], lower_bound = 0)
@variable(model, D) # incremental cycle degradation
@variable(model, C) # value of battery capacity at the end of operation
@variable(model, R) # market revenue
# binary variables
@variable(model, u2[1:T], Bin)
@variable(model, u3[1:T], Bin)

# arbitrage revenue
@constraint(model, ArbRev, R == M*sum(L.*(d-c)) )
# degradation increments
# piece-wise linear degradation opportunity value
@constraint(model, DegCost, C == M*MC*sum(d)  )
# sum of charge/discharge
@constraint(model, SumDisC[t=1:T], d[t] == d1[t] + d2[t] + d3[t] )
@constraint(model, SumChar[t=1:T], c[t] == c1[t] + c2[t] + c3[t] )
# initial SoC evolution
@constraint(model, SoCInit1, e1[1] - e01 == M*(c1[1] - d1[1]) )
@constraint(model, SoCInit2, e2[1] - e02 == M*(c2[1] - d2[1]) )
@constraint(model, SoCInit3, e1[1] - e01 == M*(c3[1] - d3[1]) )
# rest SoC evolution
@constraint(model, SoCCont1[t = 2:T], e1[t] - e1[t-1] == M*(c1[t]*eta1 - d1[t]/eta1) )
@constraint(model, SoCCont2[t = 2:T], e2[t] - e2[t-1] == M*(c2[t]*eta2 - d2[t]/eta2) )
@constraint(model, SoCCont3[t = 2:T], e3[t] - e3[t-1] == M*(c3[t]*eta3 - d3[t]/eta3) )
# final energy level
@constraint(model, Enelast, e1[T] + e2[T] +e3[T] >= ef )
# max discharge power
@constraint(model, DisMax[t=1:T], d1[t] + d2[t] + d3[t] <= P*Lp[t] )
# @constraint(model, DisMax[t=1:T], d[t] == 0*Lp[t] )
# max charge power
@constraint(model, ChrMax[t=1:T], c1[t] + c2[t] + c3[t] <= P )
# max energy level
@constraint(model, SoCMax1[t=1:T], e1[t] <= E1 )
@constraint(model, SoCMax2[t=1:T], e2[t] <= E2*u2[t] )
@constraint(model, SoCMax3[t=1:T], e3[t] <= E3*u3[t] )
# min energy level
@constraint(model, SoCMin1[t=1:T], e1[t] >= E1*u2[t] )
@constraint(model, SoCMin2[t=1:T], e2[t] >= E2*u3[t] )
@constraint(model, SoCMin3[t=1:T], e3[t] >= 0 )
# binary constraints
@constraint(model,Bin2[t=1:T], E1 - e1[t] <= 1 - u2[t] )
@constraint(model,Bin3[t=1:T], E2 - e2[t] <= 1 - u3[t] )
# warranty constraint
# @constraint(model, ETHMax, M*sum(d)/eta <= E-E*D0 )

# maximize revenue plus degradation value
@objective(model, Max, R-C)



# initialize
R_s = zeros(1, N_sim)
P_s = zeros(1, N_sim)
C_s = zeros(1, N_sim)
D_s = zeros(1, N_sim)
C_v = zeros(288, N_sim)
D_v = zeros(288, N_sim)
soc1 = zeros(288, N_sim)
soc2 = zeros(288, N_sim)
soc3 = zeros(288, N_sim)

@time begin


@printf("Optimization starts...\n")
for n = (N_sim-365):(N_sim)
# for n = (N_sim-0):(N_sim)


# update prices
local L = RTP[:,n]
# local L = collect(1:288)
local Lp = L .> 0

# update prices in constraints
for t = 1:T
    set_normalized_coefficient(ArbRev, d[t], -M*L[t] )
    set_normalized_coefficient(ArbRev, c[t], M*L[t] )
    set_normalized_rhs(DisMax[t], P*Lp[t])
    # set_normalized_rhs(DisMax[t], 0*Lp[t])
end


optimize!(model)

global R_s[n] = value(R)# objective_value(model);

global C_s[n] = value(C)
global P_s[n] = value(R)-value(C)
global D_s[n] = sum(value.(d))
global C_v[:,n] = value.(c)*M
global D_v[:,n] = value.(d)*M
global soc1[:,n] = value.(e1)
global soc2[:,n] = value.(e2)
global soc3[:,n] = value.(e3)

termination_status(model)
@printf("Finished Day %d, Cum Rev %d, Cum Profit %d, Cum Cost %d, OptStatus: %s \n", n, sum(R_s), sum(P_s), sum(C_s), termination_status(model))


end
end
# save three segment dispatch
# SoC1 = vec(reshape(soc1[:,(N_sim-364):(N_sim)],(105120,1)))
# SoC2 = vec(reshape(soc2[:,(N_sim-364):(N_sim)],(105120,1)))
# SoC3 = vec(reshape(soc3[:,(N_sim-364):(N_sim)],(105120,1)))
# df2 = DataFrame(SoC1 = SoC1, SoC2 = SoC2, SoC3 = SoC3)
# CSV.write("SoC_eta.csv", df2)
#
# # save optimal dispatch to csv
# charge = vec(reshape(C_v[:,(N_sim-364):(N_sim)],(105120,1)))
# discharge = vec(reshape(D_v[:,(N_sim-364):(N_sim)],(105120,1)))
# df = DataFrame(Charge = charge, Discharge = discharge)
# CSV.write("dispatch_PI.csv", df)
