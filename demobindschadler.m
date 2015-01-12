%% IMCORR Bindschadler Ice Stream example
%
% The IMCORR feature tracking software comes with two example image pairs.
% These are Landsat TM subscenes of a portion of Ice Stream D in West
% Antarctica (Now know as Bindschadler Ice Stream). This page shows how
% these can be tracked using the ImGRAFT toolbox.
% 


datafolder=downloadDemoData('imcorr');

A=imread(fullfile(datafolder,'conv_87.png'));
B=imread(fullfile(datafolder,'conv_89.png'));

%regular grid of points to track:

tic

[du,dv,C,Cnoise,pu,pv]=templatematch(A,B);
toc

close all
image(repmat(A,[1 1 3]),'CDataMapping','scaled') %the cdatamapping is a workaround for a bug in R2014+)
axis equal off tight ij
hold on
signal2noise=C./Cnoise;
keep=(signal2noise>2)&(C>.5);
V=(du+dv*1i)*28.5/2; %m/yr - Using imaginary values to indicate direction (for convenience).
V(~keep)=nan;
Vn=abs(V);
scatter(pu(keep),pv(keep),300,Vn(keep),'.') %colors speed
quiver(pu,pv,real(V)./Vn,imag(V)./Vn,0.2,'k') %arrows show direction. 
colormap jet
caxis([0 400])
colorbar('southoutside');
