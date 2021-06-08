load('value_dep.mat')
v_dep = v;
load('value_indep.mat')
v_idp = v;
v_dif = v_dep-v_idp;
mean(v_dif,'all')
var(v_dif,0,'all')
% mean value difference is ~10 (NYC 2016-2018)
% Assumption: indpendent cases fail to capture intertemporal dependency of day-ahead
% bias, thus result in large difference in SoC value. However, real-time
% cases doesn't have same issue, probably real-time price transition is less
% depend on date than day-ahead bias. I was wondering if there's some statistics
% to address this?
p = 20000;
hold on;
plot(v_idp(:,p))
plot(v_dep(:,p))