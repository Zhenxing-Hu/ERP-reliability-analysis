clc,close all;clear;
addpath('..\data\')
load('AEP Session1.mat');
AEP1=data(:,:,1,1,1,[133,180,350]+500);
load('AEP Session2.mat');
AEP2=data(:,:,1,1,1,[133,180,350]+500);
AEP=(AEP1+AEP2)/2;
data_temp=mean(AEP,1);
for ch=1:64
    clc;disp(ch);
    for k = 1:3
        data_temp(2,ch,1,1,1,k) = mean(bootstrp(1600,@ICC,[squeeze(AEP1(:,ch,:,:,:,k)),squeeze(AEP2(:,ch,:,:,:,k))],{'A-1'}));
    end
end

addpath('..\letswave7')
LW_init();
lwdata= FLW_load.get_lwdata('filename','AEP Session1');
lwdata.data=data_temp;
lwdata.header.name='temp';
lwdata.header.xstart=1;
lwdata.header.xstep=1;
lwdata.header.datasize=size(data_temp);
CLW_save(lwdata);
LW_init();option=[];
%% option for figure
option.inputfiles{1}='D:\paper\20210510_reliability_ASV\data&code\code\temp.lw6';
option.fig2_pos=[2017,238,700,388];

%% option.axis{1}: 133ms
option.ax{1}.name='133ms';
option.ax{1}.pos=[92,227.5288,149.3841,131.456];
option.ax{1}.style='Topograph';

%% option.axis{1}.content{1}: topo1
option.ax{1}.content{1}.name='topo1';
option.ax{1}.content{1}.type='topo';
option.ax{1}.content{1}.ep=1;
option.ax{1}.content{1}.x=[1,1];
option.ax{1}.content{1}.dim='2D';
option.ax{1}.content{1}.shrink=0.95;
option.ax{1}.content{1}.headrad=0.5;
option.ax{1}.content{1}.maplimits=[-15,15];

%% option.axis{2}: 180ms
option.ax{2}.name='180ms';
option.ax{2}.pos=[288.558,227.5288,149.3841,131.456];
option.ax{2}.style='Topograph';

%% option.axis{2}.content{1}: topo1
option.ax{2}.content{1}.name='topo1';
option.ax{2}.content{1}.type='topo';
option.ax{2}.content{1}.ep=1;
option.ax{2}.content{1}.x=[2,2];
option.ax{2}.content{1}.dim='2D';
option.ax{2}.content{1}.shrink=0.95;
option.ax{2}.content{1}.headrad=0.5;
option.ax{2}.content{1}.maplimits=[-15,15];

%% option.axis{3}: 350ms
option.ax{3}.name='350ms';
option.ax{3}.pos=[485.1159,227.5288,149.3841,131.456];
option.ax{3}.style='Topograph';

%% option.axis{3}.content{1}: topo1
option.ax{3}.content{1}.name='topo1';
option.ax{3}.content{1}.type='topo';
option.ax{3}.content{1}.ep=1;
option.ax{3}.content{1}.x=[3,3];
option.ax{3}.content{1}.dim='2D';
option.ax{3}.content{1}.shrink=0.95;
option.ax{3}.content{1}.headrad=0.5;
option.ax{3}.content{1}.maplimits=[-15,15];

%% option.axis{4}: 133ms
option.ax{4}.name='133ms';
option.ax{4}.pos=[92,43.68,149.3841,131.456];
option.ax{4}.style='Topograph';

%% option.axis{4}.content{1}: topo1
option.ax{4}.content{1}.name='topo1';
option.ax{4}.content{1}.type='topo';
option.ax{4}.content{1}.ep=2;
option.ax{4}.content{1}.x=[1,1];
option.ax{4}.content{1}.dim='2D';
option.ax{4}.content{1}.shrink=0.95;
option.ax{4}.content{1}.headrad=0.5;
option.ax{4}.content{1}.maplimits=[0,0.75];

%% option.axis{5}: 180ms
option.ax{5}.name='180ms';
option.ax{5}.pos=[288.558,43.68,149.3841,131.456];
option.ax{5}.style='Topograph';

%% option.axis{5}.content{1}: topo1
option.ax{5}.content{1}.name='topo1';
option.ax{5}.content{1}.type='topo';
option.ax{5}.content{1}.ep=2;
option.ax{5}.content{1}.x=[2,2];
option.ax{5}.content{1}.dim='2D';
option.ax{5}.content{1}.shrink=0.95;
option.ax{5}.content{1}.headrad=0.5;
option.ax{5}.content{1}.maplimits=[0,0.75];

%% option.axis{6}: 350ms
option.ax{6}.name='350ms';
option.ax{6}.pos=[485.1159,43.68,149.3841,131.456];
option.ax{6}.style='Topograph';

%% option.axis{6}.content{1}: topo1
option.ax{6}.content{1}.name='topo1';
option.ax{6}.content{1}.type='topo';
option.ax{6}.content{1}.ep=2;
option.ax{6}.content{1}.x=[3,3];
option.ax{6}.content{1}.dim='2D';
option.ax{6}.content{1}.shrink=0.95;
option.ax{6}.content{1}.headrad=0.5;
option.ax{6}.content{1}.maplimits=[0,0.75];
GLW_figure(option);
saveas(gca,'..\result\Fig4 Topograph.jpg');
delete('temp.mat');
delete('temp.lw6');
