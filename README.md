# MarkovESValuation

Analytical stochastic dynamic programming (SDP) algorithm with Markov process price model for energy storage price arbitrage, implemented in MATLAB with python for data pre-processing and Markov process model training.

## Features
* Training Markov process model for real-time price (Python) and testing real-time arbitrage with historical price data (MATLAB).
* Standard Train, Test, Markov process model for different cases in NYISO.

## Requirements

* Python 3.6+
* MATLAB
* Julia 1.6.2 and Gurobi for independent perfect information arbitrage benchmark


## Usage
Training and testing data in NYISO's four zones (NYC, LONGIL, NORTH, WEST) already include in this repo, also trained Markov process models for different cases. 

### Markov process training
Use pre_processing_main.py: 
1. Set training data set duration in line 12,13;
2. Set location in line 13,14;
3. Set model training for real-time model or DAP-RTP model in line 15;
4. Other settings are in line 37 (optional_kwargs), set summer duration, time step, price nodes gap and upper/lower bounds.

### Real-time arbitrage with real-time Markov process model
For real-time Markov model, Uuse main.m; for DAP-RTP bias Markov process model use main_DA_bias.m:
Set location, price node number and gap and cases (base, independent, with seasonal/weekly pattern), and Markov model training dataset in line 1-23.
Power-to-energy ratio (P/E), efficiency, presumed marginal discharge cost can be set in line 67-74 (main.m)/line  69-76 (main_DA_bias.m).

### Benchmark simulation
For day-ahead benchmark (BEN-DA), use main_BEN_DA.m;
For perfect information benchmark, use test_arb.jl.
