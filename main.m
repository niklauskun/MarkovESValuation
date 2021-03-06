addpath(genpath('C:\Users\wenmi\Desktop\MarkovESValuation'))

load('RTP_NYC_2010_2019.mat')
Ts = 1/12; % time step
Tp = 24/Ts; % number of timepoint
DD = 365; % select days to look back
lastDay = datetime(2019,12,31);
lambda = reshape(RTP(:,(end-DD):end),numel(RTP(:,(end-DD):end)),1); 
T = numel(lambda); % number of time steps
N = 22; % number of states
G = 10; % state gap

%% load transition matrices
pindep = 0; % price independent, 1 -> True, 0 -> False
pseason = 1; % price seasonal pattern, 1 -> True, 0 -> False
pweek = 0; % price week pattern, 1 -> True, 0 -> False
totalMatrices = 24; %total matrices number in each day
start = 2016;
stop = 2018;
location = 'NYC';

% load case matrices
if pindep == 1
    M = zeros(N,N); % initialize the transition matrices series
    for s = 1:totalMatrices % load matrices
        M(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_null_',sprintf('%d',s-1),'.csv']));
    end
else
    if pseason == 0 && pweek == 0
        M = zeros(N,N); % initialize the transition matrices series
        for s = 1:totalMatrices % load matrices
            M(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_year_',sprintf('%d',s-1),'.csv']));
        end
    elseif pseason == 1 && pweek == 0
        M1 = zeros(N,N); % initialize the summer transition matrices series
        M2 = zeros(N,N); % initialize the non-summer transition matrices series
        for s = 1:totalMatrices % load matrices
            M1(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_summer_',sprintf('%d',s-1),'.csv']));
            M2(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_nonsummer_',sprintf('%d',s-1),'.csv']));
        end
    elseif pseason == 0 && pweek == 1
        M1 = zeros(N,N); % initialize the weekday transition matrices series
        M2 = zeros(N,N); % initialize the weekend transition matrices series
        for s = 1:totalMatrices % load matrices
            M1(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_weekday_',sprintf('%d',s-1),'.csv']));
            M2(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_weekend_',sprintf('%d',s-1),'.csv']));
        end
    else
        M1 = zeros(N,N); % initialize the summer weekday transition matrices series
        M2 = zeros(N,N); % initialize the summer weekend transition matrices series
        M3 = zeros(N,N); % initialize the non-summer weekday transition matrices series
        M4 = zeros(N,N); % initialize the non-summer weekend transition matrices series
        for s = 1:totalMatrices % load matrices
            M1(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_summer_weekday_',sprintf('%d',s-1),'.csv']));
            M2(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_summer_weekend_',sprintf('%d',s-1),'.csv']));
            M3(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_nonsummer_weekday_',sprintf('%d',s-1),'.csv']));
            M4(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_','matrix_nonsummer_weekend_',sprintf('%d',s-1),'.csv']));
        end
    end
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

vEnd = zeros(Ne,1);  % generate value function samples

vEnd(1:floor(ef*100)) = 1e2; % use 100 as the penalty for final discharge level


%%
tic
v = zeros(Ne, Tp+1); % initialize the value function series
% v(1,1) is the marginal value of 0% SoC at the beginning of day 1
% V(Ne, T) is the maringal value of 100% SoC at the beginning of the last operating day
v(:,end) = vEnd; % update final value function

% initialize the value function cell for N price nodes
C = cell(1, N); 
for i = 1:N
   C{i} = v;
end
C1 = C;
C2 = C;
C3 = C;
C4 = C;

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
    if pindep == 0 && ((pseason == 1 && pweek == 0) || (pseason == 0 && pweek == 1))
        %choose transition matrix for timepoint t
        if mod(ceil(t*Ts), totalMatrices) == 0
            tM1 = M1(:,:,totalMatrices); %when mod=0, use last transition matrix
            tM2 = M2(:,:,totalMatrices); %when mod=0, use last transition matrix
        else
            tM1 = M1(:,:,mod(ceil(t*Ts),totalMatrices)); %when mode!=0
            tM2 = M2(:,:,mod(ceil(t*Ts),totalMatrices)); %when mode!=0
        end
        %calculate expected value function at timepoint t, price node i
        for i = 1:N
            viE1 = 0;
            viE2 = 0;
            for j = 1:N
                vi1 = C1{j}(:,t+1); % input value function from next timepoint at price node j
                viE1 = viE1 + tM1(i,j) * vi1; % calculate expected value function from next timepoint at price node i
                vi2 = C2{j}(:,t+1); % input value function from next timepoint at price node j
                viE2 = viE2 + tM2(i,j) * vi2; % calculate expected value function from next timepoint at price node i
            end
            lambdaNode = (i-1) *10; % calculate expected price at price node i
            vo1 = CalcValueNoUnc(lambdaNode, c, P, eta, viE1, ed, iC, iD);  % calculate value function at time point t and price node i
            vo2 = CalcValueNoUnc(lambdaNode, c, P, eta, viE2, ed, iC, iD);  % calculate value function at time point t and price node i
            C1{i}(:,t) = vo1; % record the result
            C2{i}(:,t) = vo2; % record the result
        end
    elseif pindep == 0 && pseason == 1 && pweek == 1
        %choose transition matrix for timepoint t
        if mod(ceil(t*Ts), totalMatrices) == 0
            tM1 = M1(:,:,totalMatrices); %when mod=0, use last transition matrix
            tM2 = M2(:,:,totalMatrices); %when mod=0, use last transition matrix
            tM3 = M3(:,:,totalMatrices); %when mod=0, use last transition matrix
            tM4 = M4(:,:,totalMatrices); %when mod=0, use last transition matrix
        else
            tM1 = M1(:,:,mod(ceil(t*Ts),totalMatrices)); %when mode!=0
            tM2 = M2(:,:,mod(ceil(t*Ts),totalMatrices)); %when mode!=0
            tM3 = M3(:,:,mod(ceil(t*Ts),totalMatrices)); %when mode!=0
            tM4 = M4(:,:,mod(ceil(t*Ts),totalMatrices)); %when mode!=0
        end
        %calculate expected value function at timepoint t, price node i
        for i = 1:N
            viE1 = 0;
            viE2 = 0;
            viE3 = 0;
            viE4 = 0;
            for j = 1:N
                vi1 = C1{j}(:,t+1); % input value function from next timepoint at price node j
                viE1 = viE1 + tM1(i,j) * vi1; % calculate expected value function from next timepoint at price node i
                vi2 = C2{j}(:,t+1); % input value function from next timepoint at price node j
                viE2 = viE2 + tM2(i,j) * vi2; % calculate expected value function from next timepoint at price node i
                vi3 = C3{j}(:,t+1); % input value function from next timepoint at price node j
                viE3 = viE3 + tM3(i,j) * vi3; % calculate expected value function from next timepoint at price node i
                vi4 = C4{j}(:,t+1); % input value function from next timepoint at price node j
                viE4 = viE4 + tM4(i,j) * vi4; % calculate expected value function from next timepoint at price node i
            end
            lambdaNode = (i-1) *10; % calculate expected price at price node i
            vo1 = CalcValueNoUnc(lambdaNode, c, P, eta, viE1, ed, iC, iD);  % calculate value function at time point t and price node i
            vo2 = CalcValueNoUnc(lambdaNode, c, P, eta, viE2, ed, iC, iD);  % calculate value function at time point t and price node i
            vo3 = CalcValueNoUnc(lambdaNode, c, P, eta, viE3, ed, iC, iD);  % calculate value function at time point t and price node i
            vo4 = CalcValueNoUnc(lambdaNode, c, P, eta, viE4, ed, iC, iD);  % calculate value function at time point t and price node i
            C1{i}(:,t) = vo1; % record the result
            C2{i}(:,t) = vo2; % record the result
            C3{i}(:,t) = vo3; % record the result
            C4{i}(:,t) = vo4; % record the result
        end
    else
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
end

% value function interpolation
if pindep == 0  && ((pseason == 1 && pweek == 0) || (pseason == 0 && pweek == 1))
    for i = 1:length(C)
        for t = 1:(1/Ts):Tp
            if sum(sum(C1{i}(49:end,t:t+11))) == 0
                C1{i}(:,t:t+11) = NaN;
            end
            if sum(sum(C2{i}(49:end,t:t+11))) == 0
                C2{i}(:,t:t+11) = NaN;
            end
        end
        C1{i} = fillmissing(C1{i}, 'nearest', 2);
        C2{i} = fillmissing(C2{i}, 'nearest', 2);
    end
elseif pindep == 0 && pseason == 1 && pweek == 1
    for i = 1:length(C)
        for t = 1:(1/Ts):Tp
            if sum(sum(C1{i}(49:end,t:t+11))) == 0
                C1{i}(:,t:t+11) = NaN;
            end
            if sum(sum(C2{i}(49:end,t:t+11))) == 0
                C2{i}(:,t:t+11) = NaN;
            end
            if sum(sum(C3{i}(49:end,t:t+11))) == 0
                C3{i}(:,t:t+11) = NaN;
            end
            if sum(sum(C4{i}(49:end,t:t+11))) == 0
                C4{i}(:,t:t+11) = NaN;
            end
        end
        C1{i} = fillmissing(C1{i}, 'nearest', 2);
        C2{i} = fillmissing(C2{i}, 'nearest', 2);
        C3{i} = fillmissing(C3{i}, 'nearest', 2);
        C4{i} = fillmissing(C4{i}, 'nearest', 2);
    end
else    
    for i = 1:length(C)
        for t = 1:(1/Ts):Tp
            if sum(sum(C{i}(49:end,t:t+11))) == 0
                C{i}(:,t:t+11) = NaN;
            end
        end
        C{i} = fillmissing(C{i}, 'nearest', 2);
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
 
if pindep == 0 && pseason == 1 && pweek == 0
    for t = 1:T % start from the first day and move forwards
        i = int32(1) * int32(lambda(t) < 0) + ...
            int32(22) * int32(lambda(t) >= 200) + ...        
        (idivide(lambda(t),int32(G)) + 2) * int32(lambda(t) >= 0 && lambda(t) < 200); %get price node i from lambda(t)
        tp = mod(t,int32(Tp));    
        date = lastDay - days((ceil((T-t+1)/(T/(DD+1))) - 1));
        year_start = datetime(year(date),1,1);
        day = days(date - year_start) + 1;
        vv = (124 <= day && day <= 284) * (C1{i}(:,Tp+1) * (tp == 0) + C1{i}(:,tp+1) * (tp ~= 0))...
            +(124 > day | day > 284) *  (C2{i}(:,Tp+1) * (tp == 0) + C2{i}(:,tp+1) * (tp ~= 0));
        [e, p] =  Arb_Value(lambda(t), vv, e, P, 1, eta, c, size(v,1));
        eS(t) = e; % record SoC
        pS(t) = p; % record Power
    end
elseif pindep == 0 && pseason == 0 && pweek == 1
    for t = 1:T % start from the first day and move forwards
        i = int32(1) * int32(lambda(t) < 0) + ...
            int32(22) * int32(lambda(t) >= 200) + ...        
            (idivide(lambda(t),int32(G)) + 2) * int32(lambda(t) >= 0 && lambda(t) < 200); %get price node i from lambda(t)
        tp = mod(t,int32(Tp));    
        date = lastDay - days((ceil((T-t+1)/(T/(DD+1))) - 1));
        vv = (isweekend(date) == 0) * (C1{i}(:,Tp+1) * (tp == 0) + C1{i}(:,tp+1) * (tp ~= 0))...
            +(isweekend(date) == 1) *  (C2{i}(:,Tp+1) * (tp == 0) + C2{i}(:,tp+1) * (tp ~= 0));
        [e, p] =  Arb_Value(lambda(t), vv, e, P, 1, eta, c, size(v,1));
        eS(t) = e; % record SoC
        pS(t) = p; % record Power
    end
elseif pindep == 0 && pseason == 1 && pweek == 1
    for t = 1:T % start from the first day and move forwards
        i = int32(1) * int32(lambda(t) < 0) + ...
            int32(22) * int32(lambda(t) >= 200) + ...        
            (idivide(lambda(t),int32(G)) + 2) * int32(lambda(t) >= 0 && lambda(t) < 200); %get price node i from lambda(t)
        tp = mod(t,int32(Tp));    
        date = lastDay - days((ceil((T-t+1)/(T/(DD+1))) - 1));
        year_start = datetime(year(date),1,1);
        day = days(date - year_start) + 1;
        vv = (isweekend(date) == 0 && (124 <= day && day <= 284)) * (C1{i}(:,Tp+1) * (tp == 0) + C1{i}(:,tp+1) * (tp ~= 0))...
            +(isweekend(date) == 1 && (124 <= day && day <= 284)) * (C2{i}(:,Tp+1) * (tp == 0) + C2{i}(:,tp+1) * (tp ~= 0))...
            +(isweekend(date) == 0 && (124 > day || day > 284)) * (C3{i}(:,Tp+1) * (tp == 0) + C3{i}(:,tp+1) * (tp ~= 0))...
            +(isweekend(date) == 1 && (124 > day || day > 284)) * (C4{i}(:,Tp+1) * (tp == 0) + C4{i}(:,tp+1) * (tp ~= 0));
        [e, p] =  Arb_Value(lambda(t), vv, e, P, 1, eta, c, size(v,1));
        eS(t) = e; % record SoC
        pS(t) = p; % record Power
    end
else
    for t = 1:T % start from the first day and move forwards
        i = int32(1) * int32(lambda(t) < 0) + ...
            int32(N) * int32(lambda(t) >= 200) + ...        
            (idivide(lambda(t),int32(G)) + 2) * int32(lambda(t) >= 0 && lambda(t) < 200); %get price node i from lambda(t)
        tp = mod(t,int32(Tp));    
        vv = C{i}(:,Tp+1) * (tp == 0) + C{i}(:,tp+1)* (tp ~= 0);
        [e, p] =  Arb_Value(lambda(t), vv, e, P, 1, eta, c, size(v,1));
        eS(t) = e; % record SoC
        pS(t) = p; % record Power
    end
end

ProfitOut = sum(pS.*lambda) - sum(c*pS(pS>0));
Revenue = sum(pS.*lambda);

fprintf('Profit=%e, revenue=%e',ProfitOut, Revenue)

solTimeOut = toc;
% save()

%% save figures
%if pindep == 1
%    cd('C:\Users\wenmi\Desktop\MarkovESValuation\figures\Null'); 
%elseif pindep == 0
%    cd('C:\Users\wenmi\Desktop\MarkovESValuation\figures\Year');
%end

%xs = 10;
%ys = 24;
%xblock = int32(Ne/xs);
%yblock = int32(Tp/ys);
%block = ones(xblock, yblock);

%for i = 1:N
%    Y = conv2(C{i},block,'valid');
%    Z = Y(1:xblock:end,1:yblock:end)/(double(xblock)*double(yblock));
%    h = heatmap(Z,'ColorLimits',[0 200],'Colormap',parula);
%    h.GridVisible = 'off';
%    h.Title = {'price node', num2str(i), 'SoC value'};
%    h.XLabel = 'Sizes';
%    h.YLabel = 'Colors';
%    xlabel('time')
%    ylabel('SoC')
%    saveas(h,sprintf('FIG%d%d.png',pindep,i))
%end
