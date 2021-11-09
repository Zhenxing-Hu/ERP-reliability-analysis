clc,close all;clear;
addpath('..\data\')
ERP_str={'AEP','SEP','VEP','P300'};
c1 = lines;
bar_color = [54,144,192; c1(3,:)*255;189,0,38]/255;

%% grand average
ERP_channel=[18,18,63,19];
ERP_channel_str={'Cz','Cz','Oz','Pz'};
ERP_point={[90,133,180],[110,150,245],[64,181,185],[223,345]};
ERP_color = [204,76,2;128,125,186]/255;
for ERP_idx=1:4
    load([ERP_str{ERP_idx},' Session1.mat']);
    data1=squeeze(mean(data(:,ERP_channel(ERP_idx),1,1,1,201:1500),1));
    load([ERP_str{ERP_idx},' Session2.mat']);
    data2=squeeze(mean(data(:,ERP_channel(ERP_idx),1,1,1,201:1500),1));
    figure();
    hold on;
    plot(-299:1000,data1,'color',ERP_color(1,:),'linewidth',2);
    plot(-299:1000,data2,'color',ERP_color(2,:),'linewidth',2);
    ylim = get(gca,'ylim');
    for point_idx=1:length( ERP_point{ERP_idx})
        plot([1,1]*ERP_point{ERP_idx}(point_idx),[ylim(1),ylim(2)],...
            'linestyle','--','color',bar_color(point_idx,:),'linewidth',2);
    end
    xlim([-299,1000])
    xlabel('time [ms]');
    ylabel('amplitude [\muV]');
    title([ERP_str{ERP_idx},' at Channel ',ERP_channel_str{ERP_idx}]);
    leg = legend('1st Session','2nd Session');
    saveas(gca,['..\result\Fig3 ',ERP_str{ERP_idx},' grand average.jpg'])
end