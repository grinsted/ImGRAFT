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
[pu,pv]=meshgrid(10:10:size(A,2),10:10:size(A,1));

uvA=[pu(:) pv(:)]; %these are the points we want to track.

whtemplate=10; %template size half-width
whsearch=30;   %search template half-width
super=1;
[dxy,C]=templatematch(A,B,uvA,whtemplate,whsearch,super,[0 0],true,'myncc');

close all
image(repmat(A,[1 1 3]),'CDataMapping','scaled') %the cdatamapping is a workaround for a bug in R2014+)
axis equal off tight ij
hold on
signal2noise=C(:,1)./C(:,2);
keep=(signal2noise>2)&(C(:,1)>.5);
V=dxy*28.5/2; %m/yr
Vn=sqrt(sum(V.^2,2));
scatter(uvA(keep,1),uvA(keep,2),300,Vn(keep),'.') %colors speed
quiver(uvA(keep,1),uvA(keep,2),V(keep,1)./Vn(keep),V(keep,2)./Vn(keep),0.2,'k') %arrows show direction. 
colormap jet
caxis([0 400])
colorbar('southoutside');
