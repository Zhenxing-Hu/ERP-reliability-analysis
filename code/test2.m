clc;clear;close all;
load('..\data\AEP Session1.mat');
AEP1=double(squeeze(data(:,:)));
load('..\data\AEP Session2.mat');
AEP2=double(squeeze(data(:,:)));
AEP=AEP1+AEP2;


load('..\data\SEP Session1.mat');
SEP1=double(squeeze(data(:,:)));
load('..\data\SEP Session2.mat');
SEP2=double(squeeze(data(:,:)));
SEP=SEP1+SEP2;

% r=zeros(450,450);
% for t1=1:450
%     clc;disp(t1);
%     parfor t2=1:450
%         r(t1,t2)=corr(AEP(:,t1),SEP(:,t2));
%     end
% end
% figure
% hold on;grid on;
% imagesc(r);


% unmix = pca(AEP','NumComponents',20);

[coeff1,score1,latent1]= pca(AEP'); 
[coeff2,score2,latent2]= pca(SEP');
r_train=zeros(40,100);
r_test=zeros(40,100);
for N=2:40
    clc;disp(N);
    parfor k =1:100
        [~,idx] = sort(randn(82,1));
        [A,B,r] = canoncorr(coeff1(idx(1:75),1:N),coeff2(idx(1:75),1:N));
        r_train(k,N)=r(1);
        r_test(k,N) = corr(coeff1(idx(76:82),1:N)*A(:,1),coeff2(idx(76:82),1:N)*B(:,1));
    end
end
% figure()
% hold on;
% plot(mean(r_train,1));
% plot(mean(r_test,1));
% title('score');
%% 
[coeff1,score1,latent1]= pca(AEP); 
[coeff2,score2,latent2]= pca(SEP);
r_train=zeros(40,100);
r_test=zeros(40,100);
for N=2:40
    clc;disp(N);
    for k =1:100
        [~,idx] = sort(randn(82,1));
        [A,B,r] = canoncorr(score1(idx(1:75),1:N),score2(idx(1:75),1:N));
        r_train(k,N)=r(1);
        r_test(k,N) = corr(score1(idx(76:82),1:N)*A(:,1),score2(idx(76:82),1:N)*B(:,1));
    end
end
figure()
hold on;
plot(mean(r_train,1));
plot(mean(r_test,1));
title('coef');

