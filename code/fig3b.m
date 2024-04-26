clc,close all;clear;
addpath('..\data\')
load('peak_amp_ASVP.mat');

ERP_str={'AEP','SEP','VEP','P300'};
c1 = lines;
bar_color = [54,144,192; c1(3,:)*255;189,0,38]/255;

ERP_channel=[18,18,63,19];
ERP_channel_str={'Cz','Cz','Oz','Pz'};
ERP_point={[90,133,180],[110,150,245],[64,181,185],[223,345]};
ERP_color = [204,76,2;128,125,186]/255;

for ERP_idx=1:4
    peak1 = peak_amp(:,[1,2]+(ERP_idx-1)*4);
    rel_peak(1,1) = mean(bootstrp(1600,@ICC,peak1,{'A-1'}));
    ci_peak(:,1) = bootci(1600,@ICC,peak1,{'A-1'});
    if ERP_idx<=3
        peak2 = peak_amp(:,[3,4]+(ERP_idx-1)*4);
        rel_peak(1,2) = mean(bootstrp(1600,@ICC,peak2,{'A-1'}));
        ci_peak(:,2) = bootci(1600,@ICC,peak2,{'A-1'});
    end
    
    load([ERP_str{ERP_idx},' Session1.mat']);
    data1=squeeze((data(:,ERP_channel(ERP_idx),1,1,1,501:1500)));
    load([ERP_str{ERP_idx},' Session2.mat']);
    data2=squeeze((data(:,ERP_channel(ERP_idx),1,1,1,501:1500)));
    rel_point=zeros(1,length( ERP_point{ERP_idx}),1);
    ci_point = zeros(2,length( ERP_point{ERP_idx}));
    for point_idx=1:length( ERP_point{ERP_idx})
        temp=ERP_point{ERP_idx}(point_idx);
        rel_point(1,point_idx) = mean(bootstrp(1600,@ICC,[data1(:,temp),data2(:,temp)],{'A-1'}));
        ci_point(:,point_idx) = bootci(1600,@ICC,[data1(:,temp),data2(:,temp)],{'A-1'});
    end
    if ERP_idx==1
        ERP_mean{1} = [rel_peak(1),rel_point(1),rel_point(2),rel_peak(2),rel_point(3)];
        ERP_CI{1}= [ci_peak(:,1),ci_point(:,1),ci_point(:,2),ci_peak(:,2),ci_point(:,3)];
        ERP_BarColor{1}=[1,1,3,2,2];
        ERP_xticklabels{1}={'N1','90','133','P2','180'};
    end
    if ERP_idx==2
        ERP_mean{2} = [rel_point(1),rel_peak(1),rel_point(2),rel_peak(2),rel_point(3)];
        ERP_CI{2}= [ci_point(:,1),ci_peak(:,1),ci_point(:,2),ci_peak(:,2),ci_point(:,3)];
        ERP_BarColor{2}=[3,1,1,2,2];
        ERP_xticklabels{2}={'110','N2','150','P2','245'};
    end
    if ERP_idx==3
        ERP_mean{3} = [rel_peak(1),rel_point(1),rel_point(2),rel_peak(2),rel_point(3)];
        ERP_CI{3}= [ci_peak(:,1),ci_point(:,1),ci_point(:,2),ci_peak(:,2),ci_point(:,3)];
        ERP_BarColor{3}=[1,1,3,2,2];
        ERP_xticklabels{3}={'N1','64','181','P2','185'};
    end
    if ERP_idx==4
        ERP_mean{4} = [rel_point(1),rel_peak(1),rel_point(2)];
        ERP_CI{4} = [ci_point(:,1),ci_peak(:,1),ci_point(:,2)];
        ERP_BarColor{4}=[3,1,1];
        ERP_xticklabels{4}={'223','P3','345'};
    end
    figure;hold on;
    b = bar(ERP_mean{ERP_idx});
    e=errorbar(1:length(ERP_mean{ERP_idx}), ERP_mean{ERP_idx}, ...
        ERP_mean{ERP_idx}-ERP_CI{ERP_idx}(1,:),...
        ERP_CI{ERP_idx}(2,:)-ERP_mean{ERP_idx});
    b.FaceColor = 'flat';
    for bar_idx=1:length(ERP_mean{ERP_idx})
        b.CData(bar_idx,:) = bar_color(ERP_BarColor{ERP_idx}(bar_idx),:);
    end
    e.LineStyle = 'none';
    e.Color = 'k';
    e.LineWidth = 0.5;
    set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1]);
    set(gca,'YMinorTick','on');
    set(gca,'xtick',1:length(ERP_mean{ERP_idx}));
    set(gca,'xticklabels',ERP_xticklabels{ERP_idx});
    ylabel('Reliability');
    title([ERP_str{ERP_idx},' at Channel ',ERP_channel_str{ERP_idx}]);
    saveas(gca,['..\result\Fig3 ',ERP_str{ERP_idx},' barplot.jpg'])
end