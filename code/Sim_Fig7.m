% Vary_param_trial
addpath('.\model\');
param = struct(...
      'iter_step',0.001,...
      'd',50,...  % natural frequency 
      'c',-10,... % decaying speed
      'nsub',82,...
      'C_sub_var',0.3,...
      'C_sub_mean',1,...
      'ntrial',120,... 
      'C_trial_mean',1,...
      'C_trial_var',linspace(0,2,5),...
      'jitter_sub_mean',600,...
      'jitter_sub_var',4,...
      'noise_var',500,...
      'len',1500);   
  
[ERP_ses1,ERP_ses2,grand_average] = Simulation(param, 'C_trial_var');
for idx_range = 1:size(ERP_ses1,1)
    parfor  n = 1:1500
        data = [squeeze(ERP_ses1(idx_range,n,:))...
                ,squeeze(ERP_ses2(idx_range,n,:))];
        rel(idx_range,n) = mean(bootstrp(200,@ICC,data,{'A-1'}));
        [~,var_trait(idx_range,n),var_state(idx_range,n),var_noise(idx_range,n)] = decompose_var_two_way(data);
    end
end

%%grand average
figure;hold on;box on;
linestyle = {':','--','-.','-'};

for j = 1:4
    plot(grand_average(j,201:end),'linewidth',1,'linestyle',linestyle{j},'color',[120-20*j,120-20*j,120-20*j]/255);
end
plot(grand_average(5,201:end),'-o','linewidth',1,'MarkerSize',1,'color',[120-20*j,120-20*j,120-20*j]/255);

set(gca,'xlim',[0,1300]);
set(gca,'xtick',[1,300,434,464,494,529,800,1300]);
ylim = get(gca,'ylim');
plot(gca,[434,434],[ylim(1),ylim(2)],'linestyle','--','color',[120,120,120]/255,'linewidth',1);
plot(gca,[464,464],[ylim(1),ylim(2)],'linestyle','--','color',[120,120,120]/255,'linewidth',1);
plot(gca,[494,494],[ylim(1),ylim(2)],'linestyle','--','color',[120,120,120]/255,'linewidth',1);
plot(gca,[529,529],[ylim(1),ylim(2)],'linestyle','--','color',[120,120,120]/255,'linewidth',1);
[hleg,hobj,~,~] = legend('\sigma_{trial} = 0','\sigma_{trial} = 0.5','\sigma_{trial} = 1','\sigma_{trial} = 1.5','\sigma_{trial} = 2');
set(hobj,'linewidth',2);
hobj(15).MarkerSize = 3;
ht = findobj(hobj,'type','text');
set(ht,'FontSize',12);
legend boxoff
ylim = set(gca,'ylim');
set(gca,'fontsize',12);
set(gca,'XTickLabelRotation',0);
xlabel(gca,'time [ms]');
ylabel(gca,'Amplitude [\muV]');
saveas(gca,'..\result\Fig7 ERP.jpg');

%% Reliability
figure;hold on;box on;
c(1,:)=[161,217,155]/255;
c(2,:)=[116,196,118]/255;
c(3,:)=[65,171,93]/255;
c(4,:)=[35,139,69]/255;
c(5,:)=[0,90,50]/255;
linew = [3,1.5,1.5,1,1.5];
for j = 1:5
    plot(rel(j,201:end),'linewidth',linew(j),'color',c(j,:));
end
ylabel('Reliability');
set(gca,'fontsize',9);
set(gca,'xTick',[]);
set(gca,'xlim',[0,1300]);
ylim = get(gca,'ylim');
plot([434,434],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
plot([464,464],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
plot([496,496],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
plot([529,529],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
set(gca,'ylim',ylim);
set(gca,'xtick',[1,300,800,1300]);
set(gca,'xticklabels',{'-300','0','500','1000'});
plot([0,1300],[0,0],'linewidth',1.5,'color','k');
leg=legend('\sigma_{trial} = 0','\sigma_{trial} = 0.5','\sigma_{trial} = 1','\sigma_{trial} = 1.5','\sigma_{trial} = 2');
leg.FontSize = 10;
set(gca,'ylim',[-0.3,1]);
legend boxoff;
saveas(gca,'..\result\Fig7 Reliability.jpg');

%% Var(Trait)
figure;hold on;box on;
c(1,:)=[254,237,222]/255;
c(2,:)=[253,190,133]/255;
c(3,:)=[253,141,60]/255;
c(4,:)=[230,85,13]/255;
c(5,:)=[166,54,3]/255;
linew = [3,1.5,1.5,1,1.5];
for j = 1:5
    plot(var_trait(j,201:end),'linewidth',linew(j),'color',c(j,:));
end
set(gca,'xlim',[1,1300]);
ylabel('Var(Trait)');
set(gca,'fontsize',9);
set(gca,'xTick',[]);
ylim = get(gca,'ylim');
plot([434,434],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
plot([464,464],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
plot([496,496],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
plot([529,529],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',0.2);
set(gca,'xtick',[1,300,800,1300]);
set(gca,'xticklabels',{'-300','0','500','1000'});
set(gca,'ylim',ylim);
leg=legend('\sigma_{trial} = 0','\sigma_{trial} = 0.5','\sigma_{trial} = 1','\sigma_{trial} = 1.5','\sigma_{trial} = 2');
leg.FontSize = 12;
legend boxoff;
ylabel('Var(Trait)');
set(gca,'fontsize',12);
xlabel(['time [ms]']);
saveas(gca,'..\result\Fig7 Var_trait.jpg');

%% Var(Noise)
figure;hold on;box on;
c(1,:)=[158,202,225]/255;
c(2,:)=[107,174,214]/255;
c(3,:)=[66,146,198]/255;
c(4,:)=[33,113,181]/255;
c(5,:)=[8,69,148]/255;
linew = [3,1.5,1.5,1,1.5];
for j = 1:5
    plot(var_noise(j,201:end),'linewidth',linew(j),'color',c(j,:));
end
set(gca,'xlim',[1,1300]);
ylim = get(gca,'ylim');
plot([434,434],[ylim(1),ylim(2)],'linestyle','--','color',[60,60,60]/255,'linewidth',0.3);
plot([464,464],[ylim(1),ylim(2)],'linestyle','--','color',[60,60,60]/255,'linewidth',0.3);
plot([496,496],[ylim(1),ylim(2)],'linestyle','--','color',[60,60,60]/255,'linewidth',0.3);
plot([529,529],[ylim(1),ylim(2)],'linestyle','--','color',[60,60,60]/255,'linewidth',0.3);
set(gca,'xtick',[1,300,800,1300]);
set(gca,'xticklabels',{'-300','0','500','1000'});
set(gca,'ylim',ylim);
leg=legend('\sigma_{trial} = 0','\sigma_{trial} = 0.5','\sigma_{trial} = 1','\sigma_{trial} = 1.5','\sigma_{trial} = 2');
leg.FontSize = 12;
legend boxoff;
ylabel('Var(Noise)');
set(gca,'fontsize',12);
xlabel(['time [ms]']);
saveas(gca,'..\result\Fig7 Var_noise.jpg');