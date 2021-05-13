fileList = dir('NYISO_DAP/*.csv');

Nd = numel(fileList); % number of days
Nt = 24; % data points per day
Tstep = 24*60/Nt;
Zone = 'WEST';

%%
% T = zeros(Nt, Nd);
P = zeros(288, Nd);
warning('off')
for i = 1:Nd
   [TT, PP] = NYISO_readRTP(strcat(fileList(i).folder, '/', fileList(i).name), Zone, Tstep);
   PP = reshape(transpose(repmat(PP,1,288/Nt)),[288,1]); % 288 time points in real-time
   fprintf('Reading Day %d\n', i)
   % T(:,i) = TT;
   P(:,i) = PP;
end

%%
DAP = P;
DAP(isnan(DAP)) = 0;
save('DAP_WEST_2010_2019.mat', 'DAP')