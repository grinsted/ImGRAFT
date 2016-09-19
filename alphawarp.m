function h=alphawarp(x,y,c,alpha)
%% ALPHAWARP drapes a semi-transparent and possibly distorted image over a figure.
%
% This can be used to display a variable on top of a background image.
%
% USAGE: h=alphawarp(x,y,c,alpha)
%
% EXAMPLE:
%
% P=imread('peppers.png');
% P=repmat(mean(im2double(P),3),[1 1 3]);%make a (gray) RGB-image so that it does not influence colorbar
% image(P);
% [X,Y]=meshgrid(100:300,100:300);
% X=X+Y*.1;
% Z=peaks((X-150)/100,(Y-150)/100);
% alphawarp(X,Y,Z,.5);
%
% Aslak Grinsted 2015 - http://www.glaciology.net

%extrappad is to ensure that each texture pixel is centered on the correct x,y
%value.

%todo: add input checking...
% Nans? set to typical value, and give full transparency? otherwise matlab will treat them as zeros and mess with colorbars. -Or let user deal with it.
x=extrappad(x);
y=extrappad(y);


h=surface(x,y,zeros(size(x)),c,'EdgeColor','none','FaceColor','texturemap');
if nargin<4;
    alpha=.5;
end
alpha=alpha.*(~isnan(c));
if any(alpha(:)<1)
    set(h,'FaceAlpha',  'texturemap', 'AlphaDataMapping', 'none', 'AlphaData',alpha);
end
view([0 0 1])
if nargout==0
  clear h
end


function x=extrappad(x) %ensure that pixel centers in corners are at right coordinates.

x=[ x(1,:)*1.5-x(2,:)*.5 ;0.5*(x(1:end-1,:)+x(2:end,:)) ; x(end,:)*1.5-x(end-1,:)*.5 ];
x=[ x(:,1)*1.5-x(:,2)*.5, 0.5*(x(:,1:end-1)+x(:,2:end)), x(:,end)*1.5-x(:,end-1)*.5 ];
