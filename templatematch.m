function [dxy,Cout]=templatematch(A,B,points,whtemplate,whsearch,super,dxyo,showprogress,method)
% Feature tracking by template matching
%
% [dxy,C]=templatematch(A,B,points,whtemplate,whsearch,[super,dxyo,showprogress,method])
%
% INPUTS
%    A,B: images
%    points: xy coordinates in A that should be located in B
%    whtemplate: half-[width,height] of template in A
%    whsearch: half-[width,height] of search region in B
%    super: supersampling factor (input to imresize)
%    dxyo: initial guess for x,y displacement.
%    showprogress: Boolean or cell-array of strings.
%    method: 'NCC', 'myNCC' or 'PC' (normalized cross correlation or phase correlation)
%
% Outputs:
%   dxy: displacement. [A(x,y) has moved to B(x+dx,y+dy)]
%   C: [peak-correlation,mean-abs-correlation]. (Ratio can be considered a signal to noise ratio)
%
% Methods: NCC and MyNCC are equivalent. NCC uses matlabs normxcorr2 if available
%          whereas MyNCC uses our own implementation (faster in most cases).
%
% ImGRAFT - An image georectification and feature tracking toolbox for MATLAB
% Copyright (C) 2014 Aslak Grinsted (www.glaciology.net)

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


if nargin<6
    super=1;
end
if nargin<7
    dxyo=[0 0];
end
if nargin<8
    showprogress=false;
end
if nargin<9
    method='NCC';
end

dxyo=round(dxyo);
dxyo(:,end+1)=dxyo(:,end);
if size(whtemplate,2)==1
    whtemplate(:,2)=whtemplate(:,1); %should hold [width, height.]
end
if size(whsearch,2)==1
    whsearch(:,2)=whsearch(:,1);
end
if any(whtemplate>=whsearch)
    error('Search window size must be greater than template size.')
end

Np=size(points,1);
if size(dxyo,1)==1
    dxyo(1:Np,:)=repmat(dxyo(1,:),Np,1);
end
dxy=nan(Np,2);
Cout=nan(Np,2);


switch upper(method)
    case 'NORMXCORR2'
        if exist('normxcorr2','file')>1
            matchfun=@matNCC; %TODO: use myNCC if float inputs? which is faster?
        else
            if ~isfloat(A),A=im2float(A); end
            if ~isfloat(B),B=im2float(B); end
            matchfun=@myNCC;
        end
    case {'MYNCC' 'NCC'}
        if ~isfloat(A),A=im2float(A); end
        if ~isfloat(B),B=im2float(B); end
        matchfun=@myNCC;
    case 'PC'
        matchfun=@phasecorr2;
    otherwise
        error('unknown method')
end

if islogical(showprogress)
    if ~showprogress
        showprogress=[];
    end
else
    if ~iscell(showprogress)
        if isnumeric(showprogress)
            showprogress=arrayfun(@num2str,showprogress,'uniformoutput',false);
        else
            showprogress=num2cell(showprogress);
        end
    end
end

%INITIALIZE PROGRESS FIGURE....
if ~isempty(showprogress)
    if islogical(showprogress)
        progressmsg='';
    else
        fh=figure;
        set(fh,'name','Templatematch progress','NumberTitle','off','renderer','opengl')
        hax=axes('pos',[0 0.01 0.5 0.95]);
        showimg(A);
        text(0.5,1,showprogress{1},'units','normalized','vert','bottom','fontname','courier','horiz','center')
        cc=zeros(2,size(points,1));
        hscatterA=mesh(points(:,[1 1])',points(:,[2 2])',zeros(2,size(points,1)),'mesh','column','marker','.','markersize',5,'cdata',cc); %bizarrely much faster than scatter
        colormap autumn
        caxis([0 1])
        hax(2)=axes('pos',[0.5 0.01 0.5 0.95]);
        showimg(B); hold on
        hscatterB=mesh(points(:,[1 1])'+dxyo(1),points(:,[2 2])'+dxyo(2),zeros(2,size(points,1)),'mesh','column','marker','.','markersize',5,'cdata',cc); %bizarrely much faster than scatter
        caxis([0 1])
        htext=text(1,1,'','units','normalized','vert','bottom','fontname','courier','horiz','right');
        text(0.5,1,showprogress{end},'units','normalized','vert','bottom','fontname','courier','horiz','center')
        linkaxes(hax,'xy');
        axis image, axis ij, drawnow
        hprogress=annotation(fh,'rectangle',[0 0 0 0.01],'color','none','facecolor','r');
    end
end
lastdraw=cputime;

if super==1
    resizefun=@(A,super)A;
else
    if (exist('imresize','file')>1)
        %use imresize if it is available. (requires image processing toolbox)
        resizefun=@(A,super)imresize(A,super);
    else
        if super>1
            super=2.^round(log2(super));
            resizefun=@(A,super)interp2(A,log2(super));
        else
            %         warning('image downsampling fallback for no image processing toolbox not implemented yet');
            %         resizefun=@(A,super)A;
            %         super=1;
            dwn=round(1./super);
            super=1./dwn;
            skipperfun=@(A,super)A(ceil(.5*dwn):dwn:end,floor(.5/super):dwn:end);
            %resizefun=@(A,super)skipperfun(filter2(fspecial('gaussian',[0 0]+round(2./super),1./super)),A);
            resizefun=@(A,super)skipperfun(filter2(ones(dwn)/dwn^2,A),super);
            
            %implement reshape based stacker for speed.
        end
    end
end


for ii=1:Np
    p=points(ii,:);
    
    try
        row=min(ii,size(whsearch,1));
        BB=B(p(2)+dxyo(ii,2)+(-whsearch(row,2):whsearch(row,2)),p(1)+dxyo(ii,1)+(-whsearch(row,1):whsearch(row,1)),:);
        row=min(ii,size(whtemplate,1));
        AA=A(p(2)+(-whtemplate(row,2):whtemplate(row,2)),p(1)+(-whtemplate(row,1):whtemplate(row,1)),:);
    catch
        %out of bounds
        continue
    end
    
    
    
    AA=resizefun(mean(AA,3),super); %TODO: improve edge effects of super sampling. (need 2 extra pixels for cubic) 
    BB=resizefun(mean(BB,3),super);
    
    [C,xx,yy]=matchfun(AA,BB);
    [Cmax,mix]=max(C(:));
    [mix(1),mix(2)]=ind2sub(size(C),mix);
    
    Cout(ii,2)=mean(abs(C(:))); %"noise" correlation level (we can accept that estimate even if we cannot find a good peak.)
    
    if ~(any(mix==1)||any(mix==size(C))) %do not accept any maxima on the edge of C
        %really simple/fast/crude sub pixel.  TODO: find bicubic interpolation max. (For now just super sample the imge for higher precision.)
        [xx,yy]=meshgrid(xx(mix(2)+(-1:1)),yy(mix(1)+(-1:1)));
        c=C(mix(1)+(-1:1),mix(2)+(-1:1));
        c=(c-mean(c(:)));c(c<0)=0; %best performance for landsat test images
        c=c./sum(c(:)); 
        mix(2)=sum(xx(:).*c(:));
        mix(1)=sum(yy(:).*c(:));
        
        mix=mix([2 1])./super;
        dxy(ii,:)=mix+dxyo(ii,1:2);
        Cout(ii,1)=Cmax;
    end
    
    if ~isempty(showprogress)
        if islogical(showprogress)
            if (cputime-lastdraw)>.1||(cputime<lastdraw)||(ii==Np)
                backspaces=char(zeros(size(progressmsg))+8);
                progressmsg=[uint8((1:40)<((ii/Np)*40)).*'+' ''];
                progressmsg=sprintf('Templatematch [%s]  (%5.1f %+5.1f)',progressmsg,dxy(ii,1),dxy(ii,2));
                fprintf('%s%s',backspaces,progressmsg)
                drawnow
                lastdraw=cputime;
            end
        else
            cc(:,ii)=min(max(Cout(ii,1)-Cout(ii,2),0),1);
            set(htext,'string',sprintf('%+5.1f %+5.1f ',dxy(ii,1),dxy(ii,2)))%,'units','normalized','vert','top')
            set(hprogress,'position',[0 0 ii/Np 0.01])
            set(fh,'name',sprintf('Templatematch %3.0f%%',ii*100/Np));
            if (cputime-lastdraw)>.3||(cputime<lastdraw)||(ii==Np)
                set(hscatterA,'cdata',cc);
                posB=points+dxy;
                set(hscatterB,'cdata',cc,'xdata',posB(:,[1 1])','ydata',posB(:,[2 2])');
                zlim([-1.1 .1]) %critical as otherwise matlabs clipping plane will throw out points with z=0 in older versions of matlab. BUG.
                drawnow
                lastdraw=cputime;
            end
        end
    end
    
end

if ~isempty(showprogress)
    try
        if islogical(showprogress)
            progressmsg(:)=8;
            fprintf('%s',progressmsg);
        else
            delete(htext)
            delete(hprogress)
            set(fh,'name','Templatematch: Done')
            drawnow
        end
    catch
    end
end

function [result,xx,yy]=phasecorr2(T,I)
%
%TODO-test: apply band-pass - (to remove supersampled high freqs, and
%low-passed to remove long wavelength effects (edges).
sA=size(T);sB=size(I);
sz=sA+sB-1;
myhamming=@(m)0.54 - 0.46 * cos (2 * pi * (0:m-1)' / (m-1));
ham=@(A)bsxfun(@times,bsxfun(@times,(A-mean(A(:))),myhamming(size(A,1))),myhamming(size(A,2))');
FA = conj(fft2(ham(T),sz(1),sz(2))); %2d FFT
FB = fft2(ham(I),sz(1),sz(2));
FAB = FA.*FB;
FAB = (FAB./abs(FAB));
result = real(ifft2(double(FAB)));

c=(sB-sA)/2+1; %center for zero lag
wkeep=(sB-sA)/2-3;
xx=-wkeep(2):wkeep(2);
yy=-wkeep(1):wkeep(1);
rows=mod(yy+c(1)-1,sz(1))+1;
cols=mod(xx+c(2)-1,sz(2))+1;
result=result(rows,cols);

function [C,xx,yy]=matNCC(T,B)
sT = size(T); sB = size(B);
outsize = sB + sT-1;
try
    C=normxcorr2(T,B);
catch
    C=zeros(outsize); %happens if T is constant
end
%crop to central part not affected by edges.
wkeep=(sB-sT)/2;
c=(outsize+1)/2;
C=C(c(1)+(-wkeep(1):wkeep(1)),c(2)+(-wkeep(2):wkeep(2)));
xx=-wkeep(2):wkeep(2);
yy=-wkeep(1):wkeep(1);


function [C,xx,yy]=myNCC(T,B)
%
% Alternative to NCC if no image processing toolbox. Requires single/double input.
%
sT=size(T); sB=size(B);
sz=sB+sT-1;
meanT=sum(T(:))/numel(T);
sigmaT=realsqrt(sum((T(:)-meanT).^2));
if sigmaT~=0
    fT=fft2(rot90(T,2),sz(1),sz(2));
    fB=fft2(B,sz(1),sz(2));
    C=real(ifft2(fB.*fT));
    lsumB=localsum(B,sT);
    %sigmaB=sqrt(max(localsum(B.*B,sT)-(lsumB.^2)/numel(T),0));
    %C=(C-lsumB*meanT)./(sigmaT*max(sigmaB,sigmaT/1e5)); %not 100% robust. Uses 1e5 div0 trick from D.Kroon (thanks)
    sigmaB=realsqrt(max(localsum(B.*B,sT)-(lsumB.^2)/numel(T),0)); %is the max really necessary ?
    C=(C-lsumB*meanT)./(sigmaT*sigmaB);
    C(abs(C)>1.1)=0; %this can happen if sigmaB almost 0, but we still allow C<1.1 to accomodate potential rounding issue for perfect correlation.
else
    C=zeros(sz);
end
%crop to central part not affected by edges.
wkeep=(sB-sT)/2;
c=(sz+1)/2;
C=C(c(1)+(-wkeep(1):wkeep(1)),c(2)+(-wkeep(2):wkeep(2)));
xx=-wkeep(2):wkeep(2);
yy=-wkeep(1):wkeep(1);


function lsum=localsum(A,sz)
%Fast local sum of A. Local being within the a footprint of size sz
% A = cumsum(padarray(A,sz),2);
% A = cumsum(A(:,1+sz(2):end-1)-A(:,1:end-sz(2)-1),1);
% lsum= A(1+sz(1):end-1,:)-A(1:end-sz(1)-1,:);
zp=zeros(size(A,1),sz(2)); %thanks to Matthew for spotting padarray dependency. opportunity to optimize further.
A = cumsum([zp,A,zp],2);
zp=zeros(sz(1),size(A,2));
A=[zp;A;zp];
A = cumsum(A(:,1+sz(2):end-1)-A(:,1:end-sz(2)-1),1);
lsum= A(1+sz(1):end-1,:)-A(1:end-sz(1)-1,:);


function A=im2float(A)
%im2double etc only available in image processing toolbox. this is a workaround.
if isfloat(A)
    return
end
c=class(A);
R=single([intmin(c) intmax(c)]);
A=single(A);
if R(1)~=0
    A=A-R(1);
end
A=A./diff(R);

function showimg(A) %replacement for imshow. -faster and no reliance on image processing toolbox
if isfloat(A)
    A=uint8(A*255);
else
    if strcmp(class(A),'uint16')
        A=uint8(A./256);
    end
end
% screensize=get(0,'ScreenSize');
% downsample=ceil(size(A,1)*2/screensize(3));
% if downsample>1
%     A=imresize(A,1/downsample); %TODO:remove dependency!
% end
if size(A,3)<3
    A=repmat(A,[1 1 3]);
end
[X,Y]=meshgrid([0.5 size(A,2)+0.5],[0.5 size(A,1)+0.5]);
surface(X,Y,zeros(size(X))-1,A,'EdgeColor','none','FaceColor','texturemap'); %it is much faster than using an image! (Bizarrely)
axis off tight equal image ij;
zlim([-1.1 .1]) %critical as otherwise matlabs clipping plane will throw out points with z=0 in older versions of matlab. BUG.
hold on





