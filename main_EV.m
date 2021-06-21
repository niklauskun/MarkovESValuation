addpath(genpath('C:\Users\wenmi\Desktop\MarkovESValuation'))

%% Case Setting
location = 'NYC';
startDay = datetime(2019,12,30);
endDay = datetime(2019,12,31);
yearStart = datetime(2019,1,1);
optDay = days(endDay - startDay);
numDay = days(startDay - yearStart) + 1;
Ts = 1/12; % time step
Tp = 24/Ts; % number of timepoint

% EV usage pattern
direction = 0; % 0-charge only; 1-bidirection
pluginPeriod = 1; % number of plug-in period per day

startPeriod1 = 24/Ts + 1; % start plug-in time period 1
endPeriod1 = 48/Ts; % end plug-in time period 1
eC1 = .3; % SoC consumption

startPeriod2 = 49/Ts + 1; % start plug-in time period 2
endPeriod2 = 50/Ts; % end plug-in time period 1
eC2 = .3; % SoC consumption

ed = .001; % SoC sample granularity
e0 = .0; % initial SoC 
ef = .5; % final SoC target level, use 0 if none
efp = 1e9; % use 1e9 as the penalty for final discharge level
Ne = floor(1/ed)+1; % number of SOC samples

Pr = .5; % normalized power rating wrt energy rating
P = Pr*Ts; % actual power rating taking time step size into account
eta = .9; % efficiency
c = 10; % marginal discharge cost - degradation

Nb = 12; % number of bias states (different from pre-processing, use number of states to prevent error)
Gb = 10; % bias state gap

% assert(endPeriod1 < startPeriod2,'overlap plug-in periods')
% assert(endPeriod2 - Tp < startPeriod1,'overlap plug-in periods')
%% Load Price Data
load(strcat('RTP_',location,'_2010_2019.mat'))
load(strcat('DAP_',location,'_2010_2019.mat'))
DD = 365; % select days to look back
lambda = reshape(RTP(:,(end-DD+1):end),numel(RTP(:,(end-DD+1):end)),1); 
lambda_DA = reshape(DAP(:,(end-DD+1):end),numel(DAP(:,(end-DD+1):end)),1); 
bias = lambda - lambda_DA;
T = numel(lambda); % number of time steps

%% Load Transition Matrices
totalMatrices = 24; %total matrices number in each day
start = 2016;
stop = 2018;

% load expected bias spikes
Ebs = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_expected_bias_spike.csv']));

% load case bias transition matrices
M = zeros(Nb,Nb); % initialize the transition matrices series
for s = 1:totalMatrices % load matrices
    M(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_year_',sprintf('%d',s-1),'.csv']));
end

%%
tic
qEnd = zeros(Ne,Nb,1);  % generate marginal profit function samples
qEnd(1:floor(ef/ed)+1,:) = efp; % use 1e9 as the penalty for final discharge level
q1 = zeros(Ne,Nb,endPeriod1 - startPeriod1 + 2); % initialize the marginal profit function series
q2 = zeros(Ne,Nb,endPeriod2 - startPeriod2 + 2); % initialize the marginal profit function series
q1(:,:,end) = qEnd; % update final marginal profit function
q2(:,:,end) = qEnd; % update final marginal profit function

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
h = 1:24;
H = [h,h];
%
eS = zeros(T,1); % generate the SoC series
pS = eS; % generate the power series
e = e0; % initial SoC

%% Real-time Arbitrage
for d = 1:optDay
    %% Opearation Valuation Period1
    for t = endPeriod1:-1:startPeriod1
        tp = (numDay+d-2)*Tp + t; % current time point
        tH = ceil((t)*Ts)-24; % current hour
        lambdaNode = lambda_DA(tp) + ba;
        for i = 1:Nb
            viE = (M(i,:,H(tH)) * q1(:,:,t+2-startPeriod1)')'; % calculate expected value function from next timepoint at price node i
            if direction == 0
                qo = CalcValueNoUnc_Charge(lambdaNode(i), c, P, eta, viE, ed, iC, iD);  % calculate value function at time point t and price node i
            else
                qo = CalcValueNoUnc(lambdaNode(i), c, P, eta, viE, ed, iC, iD);  % calculate value function at time point t and price node i                             
            end
            q1(:,i,t+1-startPeriod1) = qo;
        end
    end
    %% Arbitrage Period1
    for t = startPeriod1:endPeriod1
        tp = (numDay+d-2)*Tp + t; % current time point
        tH = ceil((t)*Ts)-24; % current hour
        i = int32((Nb-1)/2 + ceil(bias(tp)/Gb)); % get price node i from lambda(t)
        i = max(1,min(Nb,i));
        v(:,tp+1) = (M(i,:,H(tH)) * q1(:,:,t+2-startPeriod1)')';
        if direction == 0
            [e, p] =  Arb_Value_Charge(lambda(tp), v(:,tp+1), e, P, 1, eta, c, size(v,1));
        else
            [e, p] =  Arb_Value(lambda(tp), v(:,tp+1), e, P, 1, eta, c, size(v,1));
        end
        eS(tp) = e; % record SoC
        pS(tp) = p; % record Power
    end
        e = e - eC1; % next period initial SoC
    if pluginPeriod == 2
        %% Opearation Valuation Period2
        for t = endPeriod2:-1:startPeriod2
            tp = (numDay+d-2)*Tp + t; % current time point
            tH = ceil((t)*Ts)-24; % current hour
            lambdaNode = lambda_DA(tp) + ba;
            for i = 1:Nb
                viE = (M(i,:,H(tH)) * q2(:,:,t+2-startPeriod2)')'; % calculate expected value function from next timepoint at price node i
                if direction == 0
                    qo = CalcValueNoUnc_Charge(lambdaNode(i), c, P, eta, viE, ed, iC, iD);  % calculate value function at time point t and price node i
                else
                    qo = CalcValueNoUnc(lambdaNode(i), c, P, eta, viE, ed, iC, iD);  % calculate value function at time point t and price node i
                end
                q2(:,i,t+1-startPeriod2) = qo;
            end
        end
        %% Arbitrage Period2
        for t = startPeriod2:endPeriod2
            tp = (numDay+d-2)*Tp + t; % current time point
            tH = ceil((t)*Ts)-24; % current hour
            i = int32((Nb-1)/2 + ceil(bias(tp)/Gb)); % get price node i from lambda(t)
            i = max(1,min(Nb,i));
            v(:,tp+1) = (M(i,:,H(tH)) * q2(:,:,t+2-startPeriod2)')';
            if direction == 0
                [e, p] =  Arb_Value_Charge(lambda(tp), v(:,tp+1), e, P, 1, eta, c, size(v,1));
            else
                [e, p] =  Arb_Value(lambda(tp), v(:,tp+1), e, P, 1, eta, c, size(v,1));
            end
            eS(tp) = e; % record SoC
            pS(tp) = p; % record Power
        end
        e = e - eC2; % next period initial SoC
    end
    fprintf('Day = %d\n', d)
end

ProfitOut = sum(pS.*lambda) - sum(c*pS(pS>0));
Revenue = sum(pS.*lambda);
Discharge = sum(pS(pS>0));

fprintf('Profit = %e, Revenue = %e, Discharged = %e\n', ProfitOut, Revenue, Discharge);

solTimeOut = toc;