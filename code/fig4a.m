clc,close all;clear;
addpath('..\data\')

load('AEP Session1.mat');
AEP1=squeeze(data(:,18,1,1,1,201:1500));
load('AEP Session2.mat');
AEP2=squeeze(data(:,18,1,1,1,201:1500));
AEP=(AEP1+AEP2)/2;
t =-299:1000; 

%% Grand average waveform
[~,p,~,stats]=ttest(AEP);
color=gray(700);
figure()
hold on;box on;
for k=1:1300
    if p(k)<0.05/1000/64
        fill(gca,[0,1,1,0,0]+t(k),[-15,-15,15,15,-15]*0.99,color(floor(700-abs(stats.tstat(k))*10),:),'EdgeAlpha',0);
    end
end
plot(gca,t,mean(AEP),'linewidth',2,'color',[115,115,115]/255);
title(gca,'Grand average waveform of AEP at Cz');
set(gca,'xtick',[-300,0,133,180,350,500,1000]);
ylim = get(gca,'ylim');
plot(gca,[133,133],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot(gca,[180,180],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot(gca,[350,350],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
set(gca,'ylim',ylim);
set(gca,'xlim',[-300,1000]);
set(gca,'xticklabels',{'-300','0','133','180','350','500','1000'});
xlabel('time [ms]');
ylabel('amplitude [\muV]');
saveas(gca,'..\result\Fig4 AEP grand average.jpg')

%% test-retest reliability
rel=zeros(1,1300);
ci=zeros(2,1300);
parfor k = 1:1300
    clc;disp(k);
    rel(k) = mean(bootstrp(1600,@ICC,[AEP1(:,k),AEP2(:,k)],{'A-1'}));
    ci(:,k) = bootci(1600,@ICC,[AEP1(:,k),AEP2(:,k)],{'A-1'});
end
    
figure;hold on;box on;
plot(gca,t,rel,'color',[35,139,69]/255,'linewidth',2);
fill(gca,[t,t(end:-1:1)], [ci(1,:),ci(2,end:-1:1)], 1,'facecolor',[116,196,118]/255, 'linewidth',1.5,...
    'edgecolor', 'none', 'facealpha', 0.4);

plot(gca,[133,133],[-0.4,1],'linestyle','--','color','k','linewidth',1);
plot(gca,[180,180],[-0.4,1],'linestyle','--','color','k','linewidth',1);
plot(gca,[350,350],[-0.4,1],'linestyle','--','color','k','linewidth',1);
plot(gca,[t(1),t(end)],[0,0],'linestyle','-','linewidth',1.5,'color','k');

set(gca,'ylim',[-0.4,1])
set(gca,'ytick',[-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]);
set(gca,'xlim',[-299,1000]);
set(gca,'xtick',[-300,0,133,180,350,500,1000]);
set(gca,'xticklabels',{'-300','0','133','180','350','500','1000'});
title(gca,'Pointwise test-retest reliability of AEP');
xlabel('time [ms]');
ylabel('ICC');
saveas(gca,'..\result\Fig4 AEP reliability.jpg')

%% ICC Decomposi
Trait_AEP=zeros(1,1300);
State_AEP=zeros(1,1300);
Noise_AEP=zeros(1,1300);
parfor m = 1:1300
    M=[AEP1(:,m),AEP2(:,m)];
    [n, k] = size(M);
    SStotal = var(M(:)) *(n*k - 1);
    MSBS = var(mean(M, 2)) * k;
    MSWS = sum(var(M,0, 2)) / n;
    MSBM = var(mean(M, 1)) * n; 
    MSE = (SStotal - MSBS *(n - 1) - MSBM * (k -1))/ ((n - 1) * (k - 1));
    Trait_AEP(m) =  (MSBS - MSE)/k;
    State_AEP(m) = MSWS - MSE;
    Noise_AEP(m) = MSE;
end

figure;hold on;box on;
plot(t,Trait_AEP,'color',[241,105,19]/255,'linewidth',2);
plot(t,State_AEP,'color', [74,20,134]/255 ,'linewidth',2);
plot(t,Noise_AEP,'color',[66,146,198]/255 ,'linewidth',2);
set(gca,'xlim',[-299,1000]);
set(gca,'xtick',[-300,0,133,180,350,500,1000]);
title(gca,'Pointwise ICC Decomposition of AEP');
ylim = get(gca,'ylim');
plot(gca,[133,133],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot(gca,[180,180],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot(gca,[350,350],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
set(gca,'ylim',ylim);
ylabel('Variance');
xlabel('time [ms]');
legend('Var(Trait)','Var(State)','Var(Noise)');
set(gca,'xticklabels',{'-300','0','133','180','350','500','1000'});
saveas(gca,'..\result\Fig4 AEP ICC Decomposition.jpg')

%% Scatter plot1
figure;
hold on;box on;grid on;axis square;
plot(gca,AEP1(:,133+300),AEP2(:,133+300),'ko','MarkerFaceColor','w','MarkerSize',7);
plot(gca,AEP1(:,180+300),AEP2(:,180+300),'k*','MarkerFaceColor','w','MarkerSize',7);
I1 = ICC([AEP1(:,133+300),AEP2(:,133+300)],'A-1');
I2 = ICC([AEP1(:,180+300),AEP2(:,180+300)],'A-1');
plot(gca,[-15,40],[-15,40],'-','color','k','linewidth',1);
text(gca,0.1,0.8,['ICC_{133}=',num2str(I1,'%.2f')],'fontsize',11,'Units','normalized','color','k');
text(gca,0.1,0.9,['ICC_{180}=',num2str(I2,'%.2f')],'fontsize',11,'Units','normalized','color','k');
legend(gca,'180ms','133ms','Location','southeast');
set(gca,'xlim',[-15,40]);
set(gca,'ylim',[-15,40]);
xlabel(gca,'1st session');
ylabel(gca,'2nd session');
title(gca,'Scatter plot1');
saveas(gca,'..\result\Fig4 AEP Scatter1.jpg');

%% Scatter plot2
figure;
hold on;box on;grid on;axis square;
plot(gca,AEP1(:,350+300),AEP2(:,350+300),'ko','MarkerFaceColor','w','MarkerSize',7);
plot(gca,AEP1(:,180+300),AEP2(:,180+300),'k*','MarkerFaceColor','w','MarkerSize',7);
I1 = ICC([AEP1(:,350+300),AEP2(:,350+300)],'A-1');
I2 = ICC([AEP1(:,180+300),AEP2(:,180+300)],'A-1');
plot(gca,[-15,40],[-15,40],'-','color','k','linewidth',1);
text(gca,0.1,0.8,['ICC_{350}=',num2str(I1,'%.2f')],'fontsize',11,'Units','normalized','color','k');
text(gca,0.1,0.9,['ICC_{180}=',num2str(I2,'%.2f')],'fontsize',11,'Units','normalized','color','k');
legend(gca,'180ms','350ms','Location','southeast');
set(gca,'xlim',[-15,40]);
set(gca,'ylim',[-15,40]);
xlabel(gca,'1st session');
ylabel(gca,'2nd session');
title(gca,'Scatter plot2');
saveas(gca,'..\result\Fig4 AEP Scatter2.jpg');
