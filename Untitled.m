
% plot(eta)
% hold on
% plot(etaapp)
% xlim([0 1000])
% ylim([.6 1])
% xlabel("SoC level")
% ylabel("One-way Efficiency")
% legend('MILP linearization', 'Other - three stage efficiency')
% set(gca,'xticklabel',{'0%','20%','40%','60%','80%','100%'},'yticklabel',{'60%','70%','80%','90%','100%'},'FontSize', 16)


PI = xlsread('dispatch_PI.xlsx');
yyaxis left
% plot(-pS(1:288))
% ylim([-1/12 1/12])
% hold on
% plot(PI(1:288,1)-PI(1:288,2))
plot(lambda(1:288))
xlabel('Timepoints')
ylabel('Real-time Price ($/MWh)')
yyaxis right
plot(eS(1:288))
hold on
plot(PI(1:288,3))
ylabel('SoC Level')
ylim([0 1])
xlim([1 288])
set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'},'FontSize', 10)

legend('Real-time Price', 'DB-DEP','BEN-PF','location','south')