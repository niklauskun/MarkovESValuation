addpath(genpath('C:\Users\wenmi\Desktop\MarkovESValuation'))

location = 'NYC';
load(strcat('RTP_',location,'_2010_2019.mat'))
load(strcat('DAP_',location,'_2010_2019.mat'))
Ts = 1/12; % time step
Tp = 24/Ts; % number of timepoint
DD = 365; % select days to look back
lastDay = datetime(2019,12,31);
lambda = reshape(RTP(:,(end-DD):end),numel(RTP(:,(end-DD):end)),1); 
lambda_DA = reshape(DAP(:,(end-DD):end),numel(DAP(:,(end-DD):end)),1); 
bias = lambda - lambda_DA;
T = numel(lambda); % number of time steps
Nb = 12; % number of bias states (different from pre-processing, use number of states to prevent error)
Gb = 10; % bias state gap

lambda_n = lambda(lambda<0);
lambda_s = lambda(lambda>200); 
bias_n = bias(bias>50); 
bias_s = bias(bias<-50);
std(bias_n);
std(bias_s);

%% load transition matrices
pindep = 0; % price independent, 1 -> True, 0 -> False
pseason = 0; % price seasonal pattern, 1 -> True, 0 -> False
pweek = 0; % price week pattern, 1 -> True, 0 -> False
totalMatrices = 24; %total matrices number in each day
start = 2018;
stop = 2018;

Ebs = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_expected_bias_spike.csv']));

% load case bias transition matrices
if pindep == 1
    M = zeros(Nb,Nb); % initialize the transition matrices series
    for s = 1:totalMatrices % load matrices
        M(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_null_',sprintf('%d',s-1),'.csv']));
    end
else
    if pseason == 0 && pweek == 0
        M = zeros(Nb,Nb); % initialize the transition matrices series
        for s = 1:totalMatrices % load matrices
            M(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_year_',sprintf('%d',s-1),'.csv']));
        end
    elseif pseason == 1 && pweek == 0
        M1 = zeros(Nb,Nb); % initialize the summer transition matrices series
        M2 = zeros(Nb,Nb); % initialize the non-summer transition matrices series
        for s = 1:totalMatrices % load matrices
            M1(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_summer_',sprintf('%d',s-1),'.csv']));
            M2(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_nonsummer_',sprintf('%d',s-1),'.csv']));
        end
    elseif pseason == 0 && pweek == 1
        M1 = zeros(Nb,Nb); % initialize the weekday transition matrices series
        M2 = zeros(Nb,Nb); % initialize the weekend transition matrices series
        for s = 1:totalMatrices % load matrices
            M1(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_',sprintf('%d',Gb),'_','bias_matrix_weekday_',sprintf('%d',s-1),'.csv']));
            M2(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_',sprintf('%d',Gb),'_','bias_matrix_weekend_',sprintf('%d',s-1),'.csv']));
        end
    else
        M1 = zeros(Nb,Nb); % initialize the summer weekday transition matrices series
        M2 = zeros(Nb,Nb); % initialize the summer weekend transition matrices series
        M3 = zeros(Nb,Nb); % initialize the non-summer weekday transition matrices series
        M4 = zeros(Nb,Nb); % initialize the non-summer weekend transition matrices series
        for s = 1:totalMatrices % load matrices
            M1(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_summer_weekday_',sprintf('%d',s-1),'.csv']));
            M2(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_summer_weekend_',sprintf('%d',s-1),'.csv']));
            M3(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_nonsummer_weekday_',sprintf('%d',s-1),'.csv']));
            M4(:,:,s) = readmatrix(join([location,'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_nonsummer_weekend_',sprintf('%d',s-1),'.csv']));
        end
    end
end

% for h = 1:totalMatrices
%     zero = Nb - nnz(sum(M1(:,:,h),2));
%     fprintf('All zero transition probability at hour %d = %d \n', h, zero)
% end


% 
ba = (-50-Gb/2:Gb:50+Gb/2)';
ba(1) = Ebs(2);
ba(end) = Ebs(1);
lambda_DA_m = zeros(T,1);
for d = 1:DD
    for t = 1:Tp
        tp = (d-1)*Tp + t; % current time point
        tH = ceil((t)*Ts); % current hour
        % using previous time period bias （version A）
        if d == 1 && t == 1
            i = int32((Nb-1)/2);
        else
            i = int32((Nb-1)/2 + ceil(bias(tp-1)/Gb));
        end
        % using current time period bias （version B）
%         i = int32((Nb-1)/2 + ceil(bias(tp)/Gb));
        %%
        i = max(1,min(Nb,i));
        eb = M(i,:,tH) * ba;
        lambda_DA_m(tp) = lambda_DA(tp) + eb;
    end
end

bias_m = lambda - lambda_DA_m;

% plot(bias(1:2000))
% hold on 
% plot(bias_m(1:2000))
corrcoef(lambda_DA_m,lambda)
mean(lambda)
std(lambda)
corrcoef(lambda_DA,lambda)
lambda1 = lambda(2:end);
lambda2 = lambda(1:end-1);
corrcoef(lambda1,lambda2)

