function vis=voxelviewshed(X,Y,Z,camxyz)
% Calculate a viewshed over a DEM (fast)
%
% USAGE: vis=voxelviewshed(X,Y,Z,camxyz)
%
% INPUTS:
% X,Y,Z: input DEM (regular grid).
% camxyz: 3-element vector specifying viewpoint.
%
% OUTPUT:
%    vis: boolean visibility matrix (same size as Z)
%
% EXAMPLE:
%     [X,Y]=meshgrid(-299:300);
%     Z=abs(ifft2(ifftshift((hypot(Y,X)+1e-5).^(-2.1).*exp(rand(size(X))*2i*pi)))); %example terrain inspired by rafael.pinto
%     Z=(Z-Z(300,300))/std(Z(:)); camxyz=[0,0,0.5];
%     vis=voxelviewshed(X,Y,Z,camxyz);
%     clf;
%     surf(X,Y,Z,Z.*vis-~vis*2,'EdgeColor','none','FaceColor','interp');
%     hold on;
%     plot3(camxyz(1),camxyz(2),camxyz(3),'k*')
%     camlight
%
% v0.3 Aslak Grinsted 2014

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

if nargin==0 %run an example if no input arguments.
    warning('no input arguments. Will use test dataset.... ')
    [X,Y]=meshgrid(-299:300);
    Z=abs(ifft2(ifftshift((hypot(Y,X)+1e-5).^(-2.1).*exp(rand(size(X))*2i*pi)))); %example terrain inspired by rafael.pinto
    Z=(Z-Z(300,300))/std(Z(:)); camxyz=[0.002,0.002,0.5];
    vis=voxelviewshed(X,Y,Z,camxyz);
    clf;
    surf(X,Y,Z,Z.*vis-~vis*2,'EdgeColor','none','FaceColor','interp');
    hold on; colormap jet; camlight;
    plot3(camxyz(1),camxyz(2),camxyz(3),'k.','markersize',40)
    plot3([0 camxyz(1)],[0 camxyz(2)],[0 camxyz(3)],'k-','linewidth',3)
    clear vis
    title('monotone=hidden from view | BlackDot=camera position')
    return
end

dx=abs(X(2,2)-X(1,1));
dy=abs(Y(2,2)-Y(1,1));

sz=size(Z);
X=X(:)-camxyz(1);
Y=Y(:)-camxyz(2);
Z=Z(:)-camxyz(3);

d=sqrt(X.^2+Y.^2+Z.^2); 
x=(atan2(Y,X)+pi)/(pi*2);
y=Z./d;

[~,ix]=sortrows([round(sqrt((X/dx).^2+(Y/dy).^2)) x]); %round

loopix=find(diff(x(ix))<0);
vis=true(size(X,1),1);

maxd=max(d);%TODO: optimize
N=ceil(2*pi/(dx/maxd)); %number of points in voxel horizon

voxx=(0:N)'/N;
voxy=zeros(size(voxx))-inf;

for k=1:length(loopix)-1
    lp=ix(loopix(k)+1:loopix(k+1));
    lp=lp([end 1:end 1]);
    yy=y(lp);  xx=x(lp);
    xx(1)=xx(1)-1;xx(end)=xx(end)+1;
    vis(lp(2:end-1))=interp1q(voxx,voxy,xx(2:end-1))<yy(2:end-1);
    voxy=max(voxy,interp1q(xx,yy,voxx));
end
vis=reshape(vis,sz);

