load('RTP_NYC_2010_2019.mat')
load('DAP_NYC_2010_2019.mat')
Ts = 1/12; % time step
DD = 365; % select days to look back
lambda = reshape(RTP(:,(end-DD):end),numel(RTP(:,(end-DD):end)),1); 
lambda_DA = reshape(DAP(:,(end-DD):end),numel(DAP(:,(end-DD):end)),1); 
T = numel(lambda); % number of time steps

%%
Pr = .5; % normalized power rating wrt energy rating
P = Pr*Ts; % actual power rating taking time step size into account
eta = .9; % efficiency
c = 10; % marginal discharge cost - degradation
ed = .001; % SoC sample granularity
ef = .5; % final SoC target level, use 0 if none
Ne = floor(1/ed)+1; % number of SOC samples
e0 = .5;

vEnd = zeros(Ne,1);  % generate value function samples

vEnd(1:floor(ef*100)) = 1e2; % use 100 as the penalty for final discharge level


%%
tic
v = zeros(Ne, T+1); % initialize the value function series
% v(1,1) is the marginal value of 0% SoC at the beginning of day 1
% V(Ne, T) is the maringal value of 100% SoC at the beginning of the last operating day
v(:,end) = vEnd; % update final value function

% process index
es = (0:ed:1)';
Ne = numel(es);
% calculate soc after charge vC = (v_t(e+P*eta))
eC = es + P*eta; 
% round to the nearest sample 
iC = ceil(eC/ed)+1;
iC(iC > (Ne+1)) = Ne + 2;
iC(iC < 2) = 1;
% calculate soc after discharge vC = (v_t(e-P/eta))
eD = es - P/eta; 
% round to the nearest sample 
iD = floor(eD/ed)+1;
iD(iD > (Ne+1)) = Ne + 2;
iD(iD < 2) = 1;

for t = T:-1:1 % start from the last day and move backwards
    vi = v(:,t+1); % input value function from tomorrow
    vo = CalcValueNoUnc(lambda_DA(t), c, P, eta, vi, ed, iC, iD);
    v(:,t) = vo; % record the result 
end

tElasped = toc;

%% convert value function to 5 segments
vAvg = zeros(5,T+1);

NN = (Ne-1)/5;

for i = 1:5
   vAvg(i,:) = mean(v((i-1)*NN + (1:(NN+1)),:)); 
end


%% perform the actual arbitrage
eS = zeros(T,1); % generate the SoC series
pS = eS; % generate the power series

e = e0; % initial SoC

for t = 1:T % start from the first day and move forwards
    vv = v(:,t+1); % read the SoC value for this day
   [e, p] =  Arb_Value(lambda(t), vv, e, P, 1, eta, c, size(v,1));
   eS(t) = e; % record SoC
   pS(t) = p; % record Power
end

ProfitOut = sum(pS.*lambda) - sum(c*pS(pS>0));
Revenue = sum(pS.*lambda);
fprintf('Profit=%e, revenue=%e',ProfitOut, Revenue)
solTimeOut = toc;
