
%A primitive orthorectification example.
%
%In this approach we colorize the DEM by finding the corresponding color in
%image A. Note this approach will lead to aliasing effects when the DEM is
%coarser than the image. 

close all

%% Setup file locations and load images & data

A=imread(fullfile('data','IMG_8902.jpg'));
dem=load(fullfile('data','dem')); %load DEM
gcpA=load(fullfile('data','gcp8902.txt'));%load ground control points for image A

%% Determine camera parameters for image A
%
% Calibrate camera orientation, focal lengths and radial distortion
% Note: This is taken verbatim from the Engabreen example
% 
FocalLength=30; %mm
SensorSize=[22.0 14.8]; %mm
imgsz=size(A);
f=imgsz([2 1]).*(FocalLength./SensorSize); 
cameralocation=[446722.0 7396671.0 770.0]; 
camA=camera(cameralocation,size(A),[200 0 0]*pi/180,f); %loooking west
[camA,rmse,aic]=camA.optimizecam(gcpA(:,1:3),gcpA(:,4:5),'00000111110010000000');


%% Orthorectify by draping the photo onto the DEM using the camA projection.
%
% Note: this simple 2d interpolation will introduce some aliasing because the 
% dem is lower resolution than the photo. 
%
dem.visible=voxelviewshed(dem.X,dem.Y,dem.Z,camA.xyz);

A=im2double(A);
[uv,~,inframe]=camA.project([dem.X(:) dem.Y(:) dem.Z(:)]);
uv(~inframe|~dem.visible(:),:)=nan;

rgb=nan(size(dem.Z,1),size(dem.Z,2),1);
for jj=1:3
    rgb(:,:,jj)=reshape(interp2(A(:,:,jj),uv(:,1),uv(:,2)),size(dem.X));
end
imshow(rgb);

