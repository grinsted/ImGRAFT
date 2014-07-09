%% Batura Glacier, Karakoram
%
% CIAS is another free feature tracking software. This is a GUI written in 
% IDL by Kääb & Vollmer. Here we use ImGRAFT to track one of the CIAS example data sets using orthorectified Landsat 7 images from Batura (sample data attached below). Here's how you might track it in ImGRAFT. The single highlighted line is the one doing most of the business.
% 
% Note: this example needs the mapping toolbox in order to read the
% geo-tiffs.

datafolder=downloadDemoData('cias');
[A,Ia]=geotiffread(fullfile(datafolder,'batura_2001.tif'));
[B,Ib]=geotiffread(fullfile(datafolder,'batura_2002.tif'));
[x,y]=pixcenters(Ia,size(A));
dx=abs(x(2)-x(1));%m/pixel

%regular grid of points to track:
[pu,pv]=meshgrid(1:20:size(A,2),1:20:size(A,1));

%... but restricted to points inside this region of interest polygon
roi=[387 452;831 543;1126 899;1343 1006;1657 1022;2188 1330;...
     2437 1220;2564 1359;2483 1473;2188 1489;1693 1320;1563 1181; ...
     1061 1168;663 718;456 686;25 877;28 627;407 465];
 
ix=find(inpolygon(pu,pv,roi(:,1),roi(:,2)));
uvA=[pu(ix) pv(ix)]; %these are the points we want to track.

whtemplate=10; %template/chip size
whsearch=30; %search window size
super=1; %No super sampling 
[dxy,C]=templatematch(A,B,uvA,whtemplate,whsearch,super,[0 0],{'2001' '2002'},'myncc');

%visualize the results
close all
%turn the intensity image into an RGB image
%so that it does not interfere with colorbar:
imshow(repmat(A,[1 1 3])) 
hold on

signal2noise=C(:,1)./C(:,2);
keep=(signal2noise>2.3)&(C(:,1)>.7);
V=dxy*dx; 
Vn=sqrt(sum(V.^2,2));
scatter(uvA(keep,1),uvA(keep,2),200,Vn(keep),'.') %colors speed
quiver(uvA(keep,1),uvA(keep,2),V(keep,1)./Vn(keep),V(keep,2)./Vn(keep),0.2,'k') %arrows show direction. 
colormap jet
caxis([0 200])
hcb=colorbar('southoutside');
if ~verLessThan('matlab', '8.4.0')
    set(hcb,'limits',caxis) %workaround for a critical colorbar bug in Matlab 2014b preview... TODO: check if bug present in final release
end


