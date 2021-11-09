clc;clear;close all;
addpath('..\data\')
    
ERP_str={'AEP','SEP','VEP','P300'};
ERP_channel=[18,18,63,19];
ERP_channel_str={'Cz','Cz','Oz','Pz'};
r_std=zeros(4,1);
p_std=zeros(4,1);
r_abs_t=zeros(4,1);
p_abs_t=zeros(4,1);
for ERP_idx=1:4
    load([ERP_str{ERP_idx},' Session1.mat']);
    data1=squeeze(data(:,ERP_channel(ERP_idx),1,1,1,:));
    load([ERP_str{ERP_idx},' Session2.mat']);
    data2=squeeze(data(:,ERP_channel(ERP_idx),1,1,1,:));
    [~,p_value,~,stats]= ttest((data1+data2)/2,0); % 
    idx = find(p_value < 0.05/1500); %
    rel=zeros(1,length(data1));
    for point_idx=1:length(data1)
        rel(point_idx) = ICC([data1(:,point_idx),data2(:,point_idx)],'A-1');
    end
    
    sig_t = detrend(abs(stats.tstat(idx)))';
    sig_std = detrend(std(data1(:,idx)+data2(:,idx),1))';
    sig_rel = detrend(rel(idx))';
    
    [r_std(ERP_idx),p_std(ERP_idx)] = corr(sig_std,sig_rel,'type','spearman');
    [r_abs_t(ERP_idx),p_abs_t(ERP_idx)] = corr(sig_t,sig_rel,'type','spearman');
end
format long
disp([r_abs_t,p_abs_t,r_std,p_std]);
format short