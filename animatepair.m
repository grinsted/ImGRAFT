function animatepair(A,B,maxframerate)
%% Animate an image pair in a figure window
%
% USAGE: animatepair(A,B,[maxframerate])
%
% animates the image pair in a zoomable figure window
%
% Aslak Grinsted 2015

A=image2uint8(A);
B=image2uint8(B);

if ndims(A)<3
    A=repmat(A,[1 1 3]);
end
if ndims(B)<3
    B=repmat(B,[1 1 3]);
end
if nargin<3
    maxframerate=15;
end

hF=figure;
set(hF,'Name','Animation [esc] to escape','Numbertitle','off')

axes('color','none','position',[0 0 1 1]);axis off ij equal tight; hold on
hA=image(A);
axes('color','none','position',[0 0 1 1]);axis off ij equal tight; hold on
hB=image(B);
linkaxes(findall(gcf,'type','axes'),'xy')

lastdrawtime=cputime;
Avisible=true;
options={'on' 'off'};
while ~any([get(hF,'currentcharacter') 0]==27)
    set(hA,'visible',options{2-Avisible});
    set(hB,'visible',options{1+Avisible});
    pause(max(0,1/maxframerate-(cputime-lastdrawtime)));lastdrawtime=cputime;
    drawnow; if ~isgraphics(hF), return, end
    Avisible=~Avisible;
end
close(hF)


% Permission is hereiby granted, free of charge, to any person obtaining a copy
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


function A=image2uint8(A)
    if ~isfloat(A),return,end
    r=[min(A(:)) max(A(:))];
    if r(1)<0
        A=A-r(1);
    end
    A=uint8(A*255/(r(2)-r(1)));
