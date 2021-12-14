%% Cases Comparison
X = categorical({'RT-Idp','RT-Dep','RT-Dep-S','RT-Dep-W','DB-Idp','DB-Dep','DB-Dep-S','DB-Dep-W'});
X = reordercats(X,{'RT-Idp','RT-Dep','RT-Dep-S','RT-Dep-W','DB-Idp','DB-Dep','DB-Dep-S','DB-Dep-W'});
Y1 = [50.14 61.73 62.61 60.33 62.94 71.98 72.36 72.03];
b = bar(X,Y1,'FaceColor','flat');
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.4940 0.1840 0.5560];
b.CData(6,:) = [0.8500 0.3250 0.0980];
b.CData(7,:) = [0.9290 0.6940 0.1250];
b.CData(8,:) = [0.4940 0.1840 0.5560];
yline(100,'-','BEN-PF')
benda = yline(61.92,'--r','LineWidth',3);
benda.FontSize = 15;
set(gca,'yticklabel',{'0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%','110%'})
legend(benda,'BEN-DA')
title('NYC Profit')
saveas(gcf,'pattern.png')

%% Location/Data Size Comparison
% X = categorical({'18-18','17-18','16-18'});
% X = reordercats(X,{'18-18','17-18','16-18'});
% Y1 = [61.08 71.41 64.10;61.58 71.41 62.79;61.73 71.98 62.94];
% b = bar(X,Y1,'FaceColor','flat');
% l1 = yline(61.92,'--r','LineWidth',3);
% l2 = yline(100,'--');
% set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'},'FontSize', 20)
% legend('RT-Dep','DB-Dep','DB-Idp','BEN-DA','FontSize', 18, 'Orientation','horizontal','NumColumns',2)
% saveas(gcf,[pwd '\pictures\NYC.png'])
% 
% X = categorical({'18-18','17-18','16-18'});
% X = reordercats(X,{'18-18','17-18','16-18'});
% Y1 = [66.58 64.86 56.30;66.93 65.06 55.98;66.33 65.07 56.75];
% b = bar(X,Y1,'FaceColor','flat');
% l1 = yline(59.15,'--r','LineWidth',3);
% l2 = yline(100,'--');
% set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'},'FontSize', 20)
% legend('RT-Dep','DB-Dep','DB-Idp','BEN-DA','FontSize', 18, 'Orientation','horizontal','NumColumns',2)
% saveas(gcf,[pwd '\pictures\LONGIL.png'])
% 
% X = categorical({'18-18','17-18','16-18'});
% X = reordercats(X,{'18-18','17-18','16-18'});
% Y1 = [58.66 73.49 74.03;61.83 74.03 73.44;63.00 74.59 72.93];
% b = bar(X,Y1,'FaceColor','flat');
% l1 = yline(65.14,'--r','LineWidth',3);
% l2 = yline(100,'--');
% set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'},'FontSize', 20)
% legend('RT-Dep','DB-Dep','DB-Idp','BEN-DA','FontSize', 18, 'Orientation','horizontal','NumColumns',2)
% saveas(gcf,[pwd '\pictures\NORTH.png'])
% 
% X = categorical({'18-18','17-18','16-18'});
% X = reordercats(X,{'18-18','17-18','16-18'});
% Y1 = [75.83 78.17 73.76;75.89 78.56 73.80;75.42 77.34 73.36];
% b = bar(X,Y1,'FaceColor','flat');
% l1 = yline(72.68,'--r','LineWidth',3);
% l2 = yline(100,'--');
% set(gca,'yticklabel',{'0%','20%','40%','60%','80%','100%'},'FontSize', 20)
% legend('RT-Dep','DB-Dep','DB-Idp','BEN-DA','FontSize', 18, 'Orientation','horizontal','NumColumns',2)
% saveas(gcf,[pwd '\pictures\WEST.png'])

%% Duration/MC Comparison
% X = [0 10 30 50];
% Y = [15.3468 13.8563 10.4951 8.9827;
%     10.3162	8.4997	6.3074	5.2114;
%     6.7636	5.1422	3.5702	2.8763;
%     29.3780	21.6050	14.9360	11.7490;
%     16.9520	12.1610	8.1070	6.2430;
%     9.5890	6.6950	4.2740	3.2150;
%     ];
% p = plot(X,Y,'LineWidth',2);
% p(1).Color = [0 0.4470 0.7410];
% p(2).Color = [0.8500 0.3250 0.0980];
% p(3).Color = [0.9290 0.6940 0.1250];
% p(4).Color = [0 0.4470 0.7410];
% p(5).Color = [0.8500 0.3250 0.0980];
% p(6).Color = [0.9290 0.6940 0.1250];
% p(4).LineStyle = '--';
% p(5).LineStyle = '--';
% p(6).LineStyle = '--';
% legend('P/E = 1','P/E = 0.5','P/E = 0.25','BEN-PF P/E = 1','BEN-PF P/E = 0.5','BEN-PF P/E = 0.25','NumColumns',2)
% title('NYC Prorated Profit (1MWh Storage)')
% xlabel('Presumed Marginal Cost ($/MWh)')
% ylabel('Prorated Profit (k$)')
% set(gca,'FontSize', 14)
% saveas(gcf,'Duration MC.png')

% X = [0 10 30 50];
% Y = [12.68	59.51	123.32	174.66;
%     13.82	49.98	104.30	154.45;
%     14.50	41.92	90.35	139.42;
%     ];
% p = plot(X,Y,'LineWidth',2);
% p(1).Color = [0 0.4470 0.7410];
% p(2).Color = [0.8500 0.3250 0.0980];
% p(3).Color = [0.9290 0.6940 0.1250];
% legend('P/E = 1','P/E = 0.5','P/E = 0.25','location','northwest')
% xlabel('Presumed Marginal Cost ($/MWh)')
% ylabel('Revenue per MWh($)')
% title('NYC Revenue per MWh discharged')
% set(gca,'FontSize', 14)
% saveas(gcf,'Duration MC P.png')

%% Compare transition probabilities between training and testing data
% addpath(genpath('C:\Users\wenmi\Desktop\MarkovESValuation'))
% t = tiledlayout(2,2);
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% location = ["NYC" "LONGIL" "NORTH" "WEST"];
% for l = 1:numel(location)
% load(strcat('RTP_',location(l),'_2010_2019.mat'))
% load(strcat('DAP_',location(l),'_2010_2019.mat'))
% DD = 365; % select days to look back
% lambda = reshape(RTP(:,(end-DD+1):end),numel(RTP(:,(end-DD+1):end)),1); 
% lambda_DA = reshape(DAP(:,(end-DD+1):end),numel(DAP(:,(end-DD+1):end)),1); 
% bias = lambda - lambda_DA;
% Nb = 12;
% Gb = 10;
% totalMatrices = 24;
% 
% nps = zeros(24,1);
% nns = zeros(24,1);
% 
% biasr = reshape(bias,12,8760);
% for s = 1:totalMatrices
%     biasrs = biasr(:,s:24:end);
%     nps(s) = sum(sum(biasrs>50));
%     nns(s) = sum(sum(biasrs<-50));
% end
% 
% start = 2016;
% stop = 2018;
% ps1 = zeros(24,1);
% ns1 = zeros(24,1);
% psp1 = zeros(12,24);
% nsp1 = zeros(12,24);
% M = zeros(Nb,Nb); % initialize the transition matrices series
% for s = 1:totalMatrices % load matrices
%     M(:,:,s) = readmatrix(join([char(location(l)),'_',sprintf('%d',start),'_',sprintf('%d',stop),'_',sprintf('%d',Gb),'_','bias_matrix_year_',sprintf('%d',s-1),'.csv']));
%     ps1(s) = sum(M(:,12,s));
%     ns1(s) = sum(M(:,1,s));
%     psp1(:,s) = M(:,12,s);
%     nsp1(:,s) = M(:,1,s);
% end
% 
% ps2 = zeros(24,1);
% ns2 = zeros(24,1);
% psp2 = zeros(12,24);
% nsp2 = zeros(12,24);
% 
% Mt = zeros(Nb,Nb); % initialize the transition matrices series
% for s = 1:totalMatrices % load matrices
%     Mt(:,:,s) = readmatrix(join([char(location(l)),'_2019_2019_',sprintf('%d',Gb),'_','bias_matrix_year_',sprintf('%d',s-1),'.csv']));
%     ps2(s) = sum(Mt(:,12,s));
%     ns2(s) = sum(Mt(:,1,s));
%     psp2(:,s) = Mt(:,12,s);
%     nsp2(:,s) = Mt(:,1,s);
% end
% 
% corps = zeros(24,1);
% corns = zeros(24,1);
% 
% for s = 1:totalMatrices % load matrices
%     A = corrcoef(psp1(:,s),psp2(:,s));
%     B = corrcoef(nsp1(:,s),nsp2(:,s));
%     corps(s) = A(1,2);
%     corns(s) = B(1,2);
% end
% 
% x = 0:23;
% nexttile
% ylim([-0.5 1])
% plot(x,corps,x,corns,'LineWidth',2)
% if l == 1 || l == 3
%     set(gca,'yticklabel',{'-0.5','0','0.5','1'})
%     ylabel('Correlation Coefficients')
% else
%     set(gca,'yticklabel',[])
% end
% title(sprintf(location(l)))
% yyaxis left
% xlabel('Time (Hr)')
% hold on
% yyaxis right
% ylim([0 300])
% multibar = bar(x,[nps,nns],'EdgeColor','none','BarWidth',1);
% set(multibar(1), 'FaceColor',[0.3010, 0.7450, 0.9330])
% set(multibar(2), 'FaceColor',[0.9290 0.6940 0.1250])
% if l == 2 || l == 4
%     set(gca,'yticklabel',{'0','100','200','300'})
%     ylabel('Number of Spikes')
% else
%     set(gca,'yticklabel',[])
% end
% set(gca,'xtick', (0:4:23),'xticklabel',{'01','05','09','13','17','21'},'FontSize', 10)
% end
% 
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% lg  = legend('Positive Cor','Negative Cor','Num of Positive','Num of Negative','Orientation','Horizontal');
% lg.Layout.Tile = 'South';
% saveas(gcf,[pwd '\pictures\pcor.png'])