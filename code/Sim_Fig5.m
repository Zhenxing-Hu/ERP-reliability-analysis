addpath('.\model\');
param = struct(...
      'iter_step',0.001,...
      'd',50,...  % natural frequency 
      'c',-10,... % decaying speed
      'nsub',82,...
      'C_sub_var',0.5,...
      'C_sub_mean',1,...
      'ntrial',120,... 
      'C_trial_mean',1,...
      'C_trial_var',2,...
      'jitter_sub_mean',600,...
      'jitter_sub_var',12,...
      'noise_var',550);   
[ERP_ses1,ERP_ses2,grand_average] = Simulation(param);

parfor  t = 1:1500
    data = [ERP_ses1(t,:);...
            ERP_ses2(t,:)]';
    rel(t) = mean(bootstrp(200,@ICC,data,{'A-1'}));
    conf(t,:) = bootci(200,@ICC,data,{'A-1'});
    [~,var_trait(t),var_state(t),var_noise(t)] = decompose_var_two_way(data);
end

  
  
%%  grand average
figure;hold on;box on;
plot(grand_average(1,201:end)','linewidth',2,'color',[115,115,115]/255);
set(gca,'xtick',[]);
set(gca,'xlim',[1,1300]);
set(gca,'xtick',[1,300,434,464,529,800,1300]);
set(gca,'xticklabels',{'-300','0','134','164','229','500','1000'});
ylabel('Amplitude [\muV]');
ylim = get(gca,'ylim');
plot([434,434],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot([464,464],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot([529,529],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
set(gca,'ylim',ylim);
set(gca,'fontsize',12);
title('Grand average waveform of simulated ERP');
xlabel('time [ms]');
saveas(gca,'..\result\Fig5 grand_average_sim_ERP.jpg');

%% reliability 
figure;hold on;box on;
plot(rel(201:end),'color',[35,139,69]/255,'linewidth',2);
hold on;
x_axis = 1:1300;
x_plot =[x_axis, fliplr(x_axis)];
y_plot=[conf(201:end,1)', flipud(conf(201:end,2))'];
fill(x_plot, y_plot, 1,'facecolor', [116,196,118]/255, 'linewidth',1.5,...
    'edgecolor', 'none', 'facealpha', 0.4);
set(gca,'xlim',[1,1300]);
set(gca,'xtick',[1,300,434,464,529,800,1300]);
set(gca,'xticklabels',{'-300','0','134','164','229','500','1000'});
ylabel('Reliability');
xlabel(['time [ms]']);
title(gca,'Pointwise test-retest reliability');
set(gca,'fontsize',12);
saveas(gca,'..\result\Fig5 Rel_sim_ERP.jpg');

%% ICC decomposition
figure;hold on;box on;
c_trait = [241,105,19]/255;
c_noise = [66,146,198]/255;
c_state = [74,20,134]/255;
plot(var_trait(201:end),'color',c_trait,'linewidth',2);
plot(var_state(201:end),'color',c_state,'linewidth',2);
plot(var_noise(201:end),'color',c_noise,'linewidth',2);
set(gca,'xlim',[1,1300]);
ylim = get(gca,'ylim');
plot([434,434],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot([464,464],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
plot([529,529],[ylim(1),ylim(2)],'linestyle','--','color','k','linewidth',1);
set(gca,'ylim',ylim);
set(gca,'xtick',[1,300,434,464,529,800,1300]);
set(gca,'xticklabels',{'-300','0','134','164','229','500','1000'});
leg = legend('Var(Trait)','Var(State)','Var(Noise)');
leg.FontSize = 12;
ylabel('Variance');
title('Pointwise varaince decomposition');
set(gca,'fontsize',12);
xlabel(['time [ms]']);
saveas(gca,'..\result\Fig5 ICC_decomp_sim_ERP.jpg');