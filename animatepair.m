function animatepair(A,B)
%% Animate an image pair in a figure window
%
% USAGE: animatepair(A,B)
%

if ndims(A)<3
    A=repmat(A,[1 1 3]);
end
if ndims(B)<3
    B=repmat(B,[1 1 3]);
end

hF=figure;
set(hF,'Name','Animation [esc] to escape','Numbertitle','off')

axes('color','none','position',[0 0 1 1]);axis off ij equal tight; hold on
hA=image(A);
axes('color','none','position',[0 0 1 1]);axis off ij equal tight; hold on
hB=image(B);
linkaxes(findall(gcf,'type','axes'),'xy')

maxframerate=1/10;
lastdrawtime=cputime;
while ~any([get(hF,'currentcharacter') 0]==27)
    
    set(hA,'visible','off');
    set(hB,'visible','on');
    pause(max(0,maxframerate-(cputime-lastdrawtime)));lastdrawtime=cputime;
    drawnow;
    if ~isgraphics(hF), return, end
    
    set(hA,'visible','on');
    set(hB,'visible','off');
    pause(max(0,maxframerate-(cputime-lastdrawtime)));lastdrawtime=cputime;
    drawnow;
    if ~isgraphics(hF), return, end
end
close(hF)


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
