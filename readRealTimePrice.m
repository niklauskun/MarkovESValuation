% read all files under directory
for Year = 2019:2019
% Year = 2014;
fileList = dir(sprintf('RTP_%d/*.csv', Year));

% specify number of days
Nd = numel(fileList); % number of days
Nt = 288; % data points per day (288 for real-time price, 24 for day-ahead)
Tstep = 24*60/Nt; % step size in minutes
Zone = 'CENTRL';

%%
% T = zeros(Nt, Nd);
P = zeros(Nt, Nd); % initialize
warning('off')
for i = 1:Nd
   [TT, PP] = NYISO_readRTP(strcat(fileList(i).folder, '/', fileList(i).name), Zone, Tstep);
   fprintf('Reading Day %d\n', i)
   % T(:,i) = TT;
   P(1:numel(PP),i) = PP;
end

%%
RTP = P;
RTP(isnan(RTP)) = 0;
save(sprintf('RTP_CENTRAL_%d.mat', Year), 'RTP')

end