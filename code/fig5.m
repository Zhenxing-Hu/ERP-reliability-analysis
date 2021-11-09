clc,clear;close all;
nsub = 82;
ntrial = 80;
A=[-10*3/2,50*3/2;-50*3/2,-10*3/2];
h = 0.001;


C_sub_mean= 1;
C_sub_var = 0.1;
C_trial_mean = 1;
C_trial_var = 2;
jitter_sub_mean = 600;
jitter_sub_var = 1;
noise_var = 2000;

for_range = 1;
var_trait = zeros(for_range,1500);
var_state = zeros(for_range,1500);
var_noise = zeros(for_range,1500);
grand_average = zeros(for_range,1500);

% Pink noise
weight=ones(1999,1);
f=linspace(0,1000,1999);
f(1999:-1:1001)=f(2:1000);
weight(2:1999)=1./f(2:1999);

% System input
t=(1:2000)/1000;
a = 20;b=50;
u = a.*t.*exp(-b*t);


result=zeros(for_range,1500,nsub,ntrial);

       
jitter_var = 12;
for idx_range = 1:for_range
    jitter_sub_var = jitter_var(idx_range);
    for sub=1:nsub
        clc;disp(['nsub:',num2str(sub),'   idx_range:',num2str(idx_range)]);
        C_sub = C_sub_mean+randn(1,1)*C_sub_var;
        jitter = jitter_sub_mean+floor((rand(1,1)-0.5)*jitter_sub_var);
        for trial=1:ntrial
            result_temp=zeros(2000,2);
            C_trial = C_trial_mean + C_trial_var*randn(1,1);
            C = C_sub + C_trial;
            noise=real(ifft(fft(randn(1999,1)).*weight));
            for t=1:1999
                result_temp(t+1,:) = result_temp(t,:)+h*(result_temp(t,:)*A+noise(t)*noise_var);
                if t >jitter+500
                    result_temp(t+1,:) = result_temp(t+1,:)+ C*[u(t-jitter-500),0];%
                end
            end
            result_temp(:,1)=result_temp(:,1)-mean(result_temp(500+(1:500),1));
            result(idx_range,:,sub,trial)=result_temp(501:2000,1);
        end
    end
    exp1 =  squeeze(mean(result(idx_range,:,:,1:floor(ntrial/2)),4));
    exp2 =  squeeze(mean(result(idx_range,:,:,(1:floor(ntrial/2))+floor(ntrial/2)),4));
    grand_average(idx_range,:) = mean(squeeze(result(idx_range,:,:)),2);
    for  n = 1:1500
        data = [exp1(n,:);exp2(n,:)]';
        [~,var_trait(idx_range,n),var_state(idx_range,n),var_noise(idx_range,n)] = decompose_var_two_way(data);
    end
end

% parfor t = 1:1500
%     rel(t) = mean(bootstrp(200,@ICC,[exp1(t,:)',exp2(t,:)'],{'A-1'}));
%     conf(t,:) = bootci(200,@ICC,[exp1(t,:)',exp2(t,:)'],{'A-1'});
% end
figure;hold on;box off; axis off;
c1 = lines;
% plot(exp1(:,50),'color','g','linewidth',3);
plot(fliplr(squeeze(mean(result(1,:,1,21:50),4))),'color',[49,163,84]/255,'linewidth',3);
saveas(gca,'E:\×é»á\2021_09_23\fig\ERP_avg.tiff');

%% grand average
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