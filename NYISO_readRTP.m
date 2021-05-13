function [T_5min, Pout] = NYISO_readRTP(filename, Zone, Tstep)

%%
% pal =readtable('20200225pal.csv');
% Zone = 'N.Y.C.';
% Tstep = 5; % time step in minutes

pal = readtable(filename);

Price = table2array(pal(:,4));
ID = table2array(pal(:,2));
T = table2array(pal(:,1));

P = Price(strcmp(ID, Zone));
TT = T(strcmp(ID, Zone));
% T_5min = (TT(1):minutes(Tstep):TT(end))';
T_5min = (datetime(TT(1), 'InputFormat', 'MM/dd/yyyy HH:mm'):minutes(Tstep):datetime(TT(end), 'InputFormat', 'MM/dd/yyyy HH:mm'))';
T_5min = datetime(T_5min,'Format','MM/dd/yyyy HH:mm');

% [~,ia] = setdiff(TT, T_5min);

[~, iA, iB] = intersect(TT, T_5min);

Pout = zeros(size(T_5min));

Pout(iB) = P(iA);