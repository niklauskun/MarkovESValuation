addpath(genpath('C:\Users\wenmi\Desktop\MarkovESValuation'))

location = 'NYC';
load(strcat('RTP_',location,'_2010_2019.mat'))
load(strcat('DAP_',location,'_2010_2019.mat'))
Ts = 1/12; % time step
Tp = 24/Ts; % number of timepoint
DD = 1; % select days to look back
lastDay = datetime(2019,12,31);
lambda = reshape(RTP(:,(end-DD+1):end),numel(RTP(:,(end-DD+1):end)),1); 
lambda_DA = reshape(DAP(:,(end-DD+1):end),numel(DAP(:,(end-DD+1):end)),1); 
bias = lambda - lambda_DA;
T = numel(lambda); % number of time steps
Nb = 12; % number of bias states (different from pre-processing, use number of states to prevent error)
Gb = 10; % bias state gap

%% load transition matrices
totalMatrices = 24; %total matrices number in each day
start = 2016;
stop = 2018;

%load expected bias spikes
Ebs = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_expected_bias_spike.csv']));

% load case bias transition matrices
M = zeros(Nb,Nb); % initialize the transition matrices series
for s = 1:totalMatrices % load matrices
    M(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_year_',sprintf('%d',s-1),'.csv']));
end

%%
Pr = .5; % normalized power rating wrt energy rating
P = Pr*Ts; % actual power rating taking time step size into account
eta = .9; % efficiency
c = 10; % marginal discharge cost - degradation
ed = .001; % SoC sample granularity
ef = .5; % final SoC target level, use 0 if none
Ne = floor(1/ed)+1; % number of SOC samples
e0 = .0;

qEnd = zeros(Ne,Nb,1);  % generate value function samples

qEnd(1:floor(ef*100),:) = 1e9; % use 100 as the penalty for final discharge level

%%
tic
q = zeros(Ne,Nb,Tp+1); % initialize the value function series
% v(1,1) is the marginal value of 0% SoC at the beginning of day 1
% V(Ne, T) is the maringal value of 100% SoC at the beginning of the last operating day
q(:,:,end) = qEnd; % update final value function

% value function
v = zeros(Ne,(DD*288+1));

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
% bias index
ba = (-50-Gb/2:Gb:50+Gb/2)';
ba(1) = Ebs(2);
ba(end) = Ebs(1);

%
eS = zeros(T,1); % generate the SoC series
pS = eS; % generate the power series
e = e0; % initial SoC
lambda_t = transpose(linspace(-144,143,288));
for d = DD:-1:1
    %% valuation
    for t = Tp:-1:1
        tp = (d-1)*Tp + t; % current time point
        tH = ceil((t)*Ts); % current hour
        lambdaNode = lambda_DA(tp) + ba;
        for i = 1:Nb
            viE = (M(i,:,tH) * q(:,:,t+1)')'; % calculate expected value function from next timepoint at price node i
            qo = CalcValueNoUnc_Charge(lambdaNode(i), c, P, eta, viE, ed, iC, iD);  % calculate value function at time point t and price node i
%             maximum profit using the current bias discretization
%             ii = int32((Nb-1)/2 + ceil(bias(tp)/Gb)); % get price node i from lambda(t)
%             ii = max(1,min(Nb,ii));
%             qo = CalcValueNoUnc(lambdaNode(ii), c, P, eta, viE, ed, iC, iD);  % calculate value function at time point t and price node i
            q(:,i,t) = qo;
        end
    end
    %% abitrage
    for t = 1:Tp % start from the first time period to last time period
        tp = (d-1)*Tp + t; % current time point
        tH = ceil((t)*Ts); % current hour
        i = int32((Nb-1)/2 + ceil(bias(tp)/Gb)); % get price node i from lambda(t)
        i = max(1,min(Nb,i));
        v(:,tp) = (M(i,:,tH) * q(:,:,t+1)')';
    end
    fprintf('Day = %d\n', d)
end
v(:,end) = q(:,1,end);

%% abitrage
for d = 1:DD
    for t = 1:Tp % start from the first day and move forwards
        tp = (d-1)*Tp + t; % current time point
        [e, p] =  Arb_Value_Charge(lambda(tp), v(:,tp+1), e, P, 1, eta, c, size(v,1));
        eS(tp) = e; % record SoC
        pS(tp) = p; % record Power
    end
end

ProfitOut = sum(pS.*lambda) - sum(c*pS(pS>0));
Revenue = sum(pS.*lambda);
Discharge = sum(pS(pS>0));

fprintf('Profit = %e, Revenue = %e, Discharged = %e\n', ProfitOut, Revenue, Discharge);

solTimeOut = toc;