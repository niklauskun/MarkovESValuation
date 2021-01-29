cd 'C:\Users\wenmi\Desktop\MarkovESValuation'
load('RTP_NYC_2010_2019.mat')
Ts = 1/12; % time step
Tp = 24/Ts; % number of timepoint
DD = 365; % select days to look back
lambda = reshape(RTP(:,(end-DD):end),numel(RTP(:,(end-DD):end)),1); 
T = numel(lambda); % number of time steps
N = 22; %number of states
G = 10; %state gap

%% load transition matrices
oldFolder = cd('C:\Users\wenmi\Desktop\MarkovESValuation\transition_matrix'); %transition matrix folder (will add argument for diff scenario)
matrices = dir('*');
matrices(1:2) = [];
totalMatrices = numel(matrices); %total matrices number

M = zeros(N,N); % initialize the transition matrices series

for s = 1:totalMatrices % load all matrices
    M(:,:,s) = readmatrix(join(['matrix',sprintf('%d',s-1),'.csv']));
end

cd(oldFolder)

%%
Pr = .5; % normalized power rating wrt energy rating
P = Pr*Ts; % actual power rating taking time step size into account
eta = .9; % efficiency
c = 10; % marginal discharge cost - degradation
ed = .001; % SoC sample granularity
ef = .0; % final SoC target level, use 0 if none
Ne = floor(1/ed)+1; % number of SOC samples
e0 = .5;

vEnd = zeros(Ne,1);  % generate value function samples

vEnd(1:floor(ef*100)) = 1e2; % use 100 as the penalty for final discharge level


%%
tic
v = zeros(Ne, Tp+1); % initialize the value function series
% v(1,1) is the marginal value of 0% SoC at the beginning of day 1
% V(Ne, T) is the maringal value of 100% SoC at the beginning of the last operating day
v(:,end) = vEnd; % update final value function

C = cell(1, N); % initialize the value function cell for N price nodes
for i = 1:N
   C{i} = v;
end

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

for t = Tp:-1:1 % start from the last timepoint and move backwards
    %choose transition matrix for timepoint t
    if mod(ceil(t*Ts), totalMatrices) == 0
        tM = M(:,:,totalMatrices); %when mod=0, use last transition matrix
    else
        tM = M(:,:,mod(ceil(t*Ts),totalMatrices)); %when mode!=0
    end
    %calculate expected value function at timepoint t, price node i
    for i = 1:N
        viE = 0;
        for j = 1:N
            vi = C{j}(:,t+1); % input value function from next timepoint at price node j
            viE = viE + tM(i,j) * vi; % calculate expected value function from next timepoint at price node i
        end
        lambdaNode = (i-1) *10; % calculate expected price at price node i
        vo = CalcValueNoUnc(lambdaNode, c, P, eta, viE, ed, iC, iD);  % calculate value function at time point t and price node i
        C{i}(:,t) = vo; % record the result 
    end
end

tElasped = toc;

%% convert value function to 5 segments
%vAvg = zeros(5,T+1);

%NN = (Ne-1)/5;

%for i = 1:5
%   vAvg(i,:) = mean(v((i-1)*NN + (1:(NN+1)),:)); 
%end


%% perform the actual arbitrage
eS = zeros(T,1); % generate the SoC series
pS = eS; % generate the power series

e = e0; % initial SoC

for t = 1:T % start from the first day and move forwards
    i = idivide(lambda(t),int16(G)) + 2; %get price node i from lambda(t)
    if lambda(t) < 0
        i = int16(1);
    elseif lambda(t) >= 200
        i = int16(22);
    end
    tp = mod(t,int16(Tp));
    if tp == 0
        vv = C{i}(:,Tp+1);
    else
        vv = C{i}(:,tp+1); % read the SoC value for this day
    end
   [e, p] =  Arb_Value(lambda(t), vv, e, P, 1, eta, c, size(v,1));
   eS(t) = e; % record SoC
   pS(t) = p; % record Power
end

ProfitOut = sum(pS.*lambda) - sum(c*pS(pS>0));
solTimeOut = toc;
