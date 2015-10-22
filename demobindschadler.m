%% Bindschadler Ice Stream example
%
% The IMCORR feature tracking software comes with two example image pairs.
% These are Landsat TM subscenes of a portion of Ice Stream D in West
% Antarctica (Now know as Bindschadler Ice Stream). This page shows how
% these can be tracked using the ImGRAFT toolbox.
% 

datafolder=downloadDemoData('imcorr');

A=imread(fullfile(datafolder,'conv_87.png'));
B=imread(fullfile(datafolder,'conv_89.png'));

tic
[du,dv,C,Cnoise,pu,pv]=templatematch(A,B);
toc

close all
image(repmat(A,[1 1 3]),'CDataMapping','scaled') %the cdatamapping is a workaround for a bug in R2014+)
axis equal off tight ij
hold on
signal2noise=C./Cnoise;
keep=(signal2noise>3.5);
V=(du+dv*1i)*28.5/2; %m/yr - Using imaginary values to indicate direction (for convenience).
Vn=abs(V);
alphawarp(pu,pv,Vn,.2+keep*.5)
quiver(pu(keep),pv(keep),real(V(keep))./Vn(keep),imag(V(keep))./Vn(keep),0.2,'k') %arrows show direction. 
caxis([0 400])
colorbar('eastoutside');
