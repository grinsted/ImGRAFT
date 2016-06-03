function h=showimg(x,y,A)
% Display an image
%
% Usage h=showimg(x,y,A)
%
%
% Aslak Grinsted 2016

if ~isa(A,'uint8')
    A=A-min(A(:));
    A=uint8(A*255.49./max(A(:)));
end
if size(A,3)==1
    A=repmat(A,[1 1 3]); %this is to ensure that it does not interfere with CAXIS.
end

h=image(A,'XData',x,'YData',y,'CDataMapping','scaled'); %the cdatamapping is a workaround for a bug in R2014+
axis equal tight xy off
if nargout==0
    clearvars h
end
