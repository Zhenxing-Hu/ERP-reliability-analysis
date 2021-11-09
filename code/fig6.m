clc,clear;close all;
nsub = 82;
ntrial = 120;
A=[-10,50;-50,-10];
h = 0.001;


C_sub_mean= 1;
C_sub_var = 0.5;
C_trial_mean = 1;
C_trial_var = 2;
jitter_sub_mean = 600;
jitter_sub_var = 1;
noise_var = 500;

for_range = 5;
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

jitter_var = [0,1,2,10,20];
for idx_range = 1:for_range
    jitter_sub_var = jitter_var(idx_range);
    for sub=1:nsub
        clc;disp(['nsub:',num2str(sub),'   idx_range:',num2str(idx_range)]);
        C_sub = C_sub_mean+randn(1,1)*C_sub_var;
        jitter = jitter_sub_mean+floor((rand(1,1)-0.5)*2*jitter_sub_var);
        parfor trial=1:ntrial
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
    parfor  n = 1:1500
        data = [exp1(n,:);exp2(n,:)]';
        [~,var_trait(idx_range,n),var_state(idx_range,n),var_noise(idx_range,n)] = decompose_var_two_way(data);
    end
end

rel = var_trait./(var_trait + var_state + var_noise);

%% Grand average
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
[hleg,hobj,~,~] = legend('\tau_{s} = 0','\tau_{s} = 1','\tau_{s} = 2','\tau_{s} = 10','\tau_{s} = 20');

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
saveas(gca,'..\result\Fig6 ERP.jpg');

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
leg=legend('\tau_{s} = 0','\tau_{s} = 1','\tau_{s} = 2','\tau_{s} = 10','\tau_{s} = 20');
leg.FontSize = 10;
set(gca,'ylim',[-0.3,1]);
legend boxoff;
saveas(gca,'..\result\Fig6 Reliability.jpg');

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
leg=legend('\tau_{s} = 0','\tau_{s} = 1','\tau_{s} = 2','\tau_{s} = 10','\tau_{s} = 20');
leg.FontSize = 12;
legend boxoff;
ylabel('Var(Trait)');
set(gca,'fontsize',12);
xlabel(['time [ms]']);
saveas(gca,'..\result\Fig6 Var_trait.jpg');

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
leg=legend('\tau_{s} = 0','\tau_{s} = 1','\tau_{s} = 2','\tau_{s} = 10','\tau_{s} = 20');
leg.FontSize = 12;
legend boxoff;
ylabel('Var(Noise)');
set(gca,'fontsize',12);
xlabel(['time [ms]']);
saveas(gca,'..\result\Fig6 Var_noise.jpg');


