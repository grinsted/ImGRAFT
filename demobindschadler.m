
datafolder=downloadDemoData('imcorr');

A=imread(fullfile(datafolder,'conv_87.png'));
B=imread(fullfile(datafolder,'conv_89.png'));

%regular grid of points to track:
[pu,pv]=meshgrid(10:20:size(A,2),10:20:size(A,1));

uvA=[pu(:) pv(:)]; %these are the points we want to track.

whtemplate=10; %template size half-width
whsearch=30;   %search template half-width
super=1;
[dxy,C]=templatematch(A,B,uvA,whtemplate,whsearch,super,[0 0],{'1987' '1989'},'myncc');

close all
imshow(repmat(A,[1 1 3]))
hold on
signal2noise=C(:,1)./C(:,2);
keep=(signal2noise>2)&(C(:,1)>.5);
V=dxy*28.5/2; %m/yr
Vn=sqrt(sum(V.^2,2));
scatter(uvA(keep,1),uvA(keep,2),300,Vn(keep),'.') %colors speed
quiver(uvA(keep,1),uvA(keep,2),V(keep,1)./Vn(keep),V(keep,2)./Vn(keep),0.2,'k') %arrows show direction. 
colormap jet
caxis([0 400])
hcb=colorbar('southoutside');
if ~verLessThan('matlab', '8.4.0')
    set(hcb,'limits',caxis) %workaround for a critical colorbar bug in Matlab 2014b preview... TODO: check if bug present in final release
end

