close all

Int1 = imread('int1.jpg');
Int2 = imread('int2.jpg');


X1 = double(im2gray(Int1));
X2 = double(im2gray(Int2));

nx1 = size(X1,1); ny1 = size(X1,2);
nx2 = size(X2,1); ny2 = size(X2,2);


%%%%%%%%%%%%%%%% SVD ANALYSIS %%%%%%%%%%%%%%%%%%%%%%

[U1,S1,V1] = svd(X1);
[U2,S2,V2] = svd(X2);

sv1 = diag(S1);
sv2 = diag(S2);
cs1 = cumsum(sv1)/sum(sv1);
cs2 = cumsum(sv2)/sum(sv2);

f1=figure;
set(f1,"Position",[ 366   689   995   454])
subplot(1,2,1),semilogy(1:50,sv1(1:50),'.-',1:50,sv2(1:50),'.-')
legend('No waves','With waves')
title('S-matrix values')
subplot(1,2,2), plot(1:nx1,cs1, 1:nx2,cs2)
title('Cumulative S-matrix values')

modes = 3;
Xapprox1 = U1(:,1:modes)*S1(1:modes,1:modes)*V1(:,1:modes)';
Xapprox2 = U2(:,1:modes)*S2(1:modes,1:modes)*V2(:,1:modes)';

f=figure;

tit = sprintf('Number of SVD modes: %d', modes);
set(f,"Position",[ 133         849        1594         420])
colormap gray

t=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile, imagesc(X1), title('original'), axis off;%, axis equal;
nexttile, imagesc(X2), title('original'),axis off;%, axis equal;
nexttile, imagesc(Xapprox1), title(tit),axis off;%, axis equal;
nexttile, imagesc(Xapprox2), title(tit),axis off;%, axis equal;

fh1 = figure;
bins = 10;
imhist(Int1,bins); title('No waves -- greyscale histogram')
fh2 = figure;
imhist(Int2,bins);title('With waves -- greyscale histogram')

%%%%%%%%%%% CROSS CORRELATION ANALYSIS %%%%%%%%%%%%%%%%
% maxlag = 20;
% step=50;
% corrs1 = zeros(ny1-step,2*maxlag+1);
% corrs2 = zeros(ny2-step,2*maxlag+1);
% 
% 
% for i=1:ny1-step
%     corrs1(i,:) = xcorr(X1(:,i),X1(:,i+step),maxlag,'normalized');
%     corrs2(i,:) = xcorr(X2(:,i),X2(:,i+step),maxlag,'normalized');
% end
% 
% xcavg1 = mean(corrs1,1);
% xcavg2 = mean(corrs2,1);
% 
% plot(-maxlag:maxlag,xcavg1,'-k.');
% hold on
% plot(-maxlag:maxlag,xcavg2,'-r.');
% hold off
% 
% 
% %zero-lag correlation falloff
% maxstep = 100;
% clag01=zeros(1,maxstep);  clag02=zeros(1,maxstep);
% 
% for s = 1:maxstep
%   for i = ny1-maxstep
%       clag01(s) = xcorr(X1(:,i),X1(:,i+s),0,'normalized');
%       clag02(s) = xcorr(X2(:,i),X2(:,i+s),0,'normalized');
%   end
% end
% 
% close all;
% plot(1:maxstep,clag01,1:maxstep,clag02)


%%%%%%%%%%%%%%  IMAGE HISTOGRAMS %%%%%%%%%%%

% tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
% nexttile, imhist(Int1);
% nexttile, imhist(Int2);



