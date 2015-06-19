function [du, dv, peakCorr, meanAbsCorr,pu,pv]=templatematch(A,B,varargin)
%% Feature tracking by template matching
%
% 
%
% SYNTAX:
%    [du, dv, peakCorr, meanAbsCorr, pu, pv] = templatematch(A,B[,pu,pv][,parameter-value-pairs])
%
% INPUTS
%    A,B: images
%    pu,pv: pixel coordinates in A that should be located in B. (Default is a regular grid)
%
%
% NAMED PARAMETERS: 
%    TemplateWidth,TemplateHeight: Size of templates in A (Default: 21). 
%    SearchWidth,SearchHeight: Size of search region in B (Default: TemplateWidth+40).
%    SuperSample: super sampling factor of input images for improved subpixel accuracy. (default=1)
%    Initialdu,Initialdv: initial guess of the displacement between A & B
%    Super: supersampling factor (input to imresize)
%    ShowProgress: Boolean or cell-array of strings.
%                  true (default) is used for a text progress bar.
%                  A cell of strings is used to name the A & B images in a progress figure.
%    Method: 'NCC'(default), 'NORMXCORR2' or 'PC' (normalized cross correlation or phase correlation)
%
%  * The template/search height is assumed to be the same as corresponding
%  widths if not explicitly specified. 
%
% OUTPUTS:
%   du,dv: displacement of each point in pu,pv. [A(pu,pv) has moved to B(pu+du,pv+dv)]
%   peakCorr: correlation coefficient of the matched template. 
%   meanAbsCorr: The mean absolute correlation coefficitent over the search
%                region is an estimate of the noise level.
%   pu,pv: actual pixel centers of templates in A may differ from inputs because of rounding. 
%
%
% ImGRAFT - An image georectification and feature tracking toolbox for MATLAB
% Copvright (C) 2014 Aslak Grinsted (www.glaciology.net)

% Permission is hereby granted, free of charge, to any person obtaining a copv
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copv, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copvright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

%templatematch(A,B,pu,pv,whtemplate,whsearch,super,duo,dvo,ShowProgress,Method)
%              1 2 3  4   5            6       7   8    9   10           11

p=inputParser;
p.FunctionName='templatematch';
try %older versions of matlab did not support partialmatching.
    p.CaseSensitive=false;
    p.StructExpand=true;
    p.PartialMatching=true;
catch
end
p.addOptional('pu',[],@isnumeric);
p.addOptional('pv',[],@isnumeric);
p.addParameter('TemplateWidth',41,@isnumeric); 
p.addParameter('TemplateHeight',[],@isnumeric);
p.addParameter('SearchWidth',[],@isnumeric);
p.addParameter('SearchHeight',[],@isnumeric);
p.addParameter('Initialdu',0,@isnumeric);
p.addParameter('Initialdv',0,@isnumeric);
p.addParameter('ShowProgress',true); 
p.addParameter('Method','NCC',@ischar);
p.addParameter('SuperSample',1);

p.parse(varargin{:});
R=p.Results;

if all(size(R.pu))~=all(size(R.pv))
    error('imgraft:inputerror', 'pu and pv must be same size.');
end

if isempty(R.TemplateHeight), R.TemplateHeight=R.TemplateWidth; end;
if isempty(R.SearchWidth), R.SearchWidth = R.TemplateWidth+60; end;
if isempty(R.SearchHeight), R.SearchHeight = R.SearchWidth; end;

if isempty(R.pu)
    dpu=max(R.TemplateWidth(1)/2,size(A,2)/100); %dont auto generate excessive number of points (idea from: Chad Greene)
    dpv=max(R.TemplateHeight(1)/2,size(A,1)/100); 
    R.pu=R.SearchWidth(1)/2 : dpu : size(A,2)-R.SearchWidth(1)/2;
    R.pv=R.SearchHeight(1)/2: dpv : size(A,1)-R.SearchHeight(1)/2;
    [R.pu,R.pv]=meshgrid(R.pu,R.pv);
end

%TODO: make better dimension checking!
if ((numel(R.Initialdu)>1) || (numel(R.Initialdv)>1)) && ((numel(R.Initialdu)~=numel(R.pu)) || (numel(R.Initialdv)~=numel(R.pv))) 
    error('imgraft:inputerror', 'size of Initialdu and Initialdv must match pu and pv');
end
if ((numel(R.TemplateWidth)>1) || (numel(R.TemplateHeight)>1)) && ((numel(R.TemplateWidth)~=numel(R.pu)) || (numel(R.TemplateHeight)~=numel(R.pv))) 
    error('imgraft:inputerror', 'number of elements in TemplateWidth/height must match pu and pv');
end
if ((numel(R.SearchWidth)>1) || (numel(R.SearchHeight)>1)) && ((numel(R.SearchWidth)~=numel(R.pu)) || (numel(R.SearchHeight)~=numel(R.pv))) 
    error('imgraft:inputerror', 'number of elements in SearchWidth/height must match pu and pv');
end


R.Initialdu=round(R.Initialdu);
R.Initialdv=round(R.Initialdv);

if any(R.TemplateWidth>=R.SearchWidth)||any(R.TemplateHeight>=R.SearchHeight)
    error('imgraft:inputerror','Search window size must be greater than template size.')
end

Np=numel(R.pu);


du=nan(size(R.pu));
dv=nan(size(R.pu));
peakCorr=nan(size(R.pu));
meanAbsCorr=nan(size(R.pu));
pu=R.pu;
pv=R.pv;


switch upper(R.Method)
    case 'NORMXCORR2'
        if license('test','image_toolbox') %exist('normxcorr2','file')>1
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
        error('imgraft:inputerror','unknown Method: %s',R.Method)
end

if islogical(R.ShowProgress)
    if ~R.ShowProgress
        R.ShowProgress=[];
    end
else
    if ~iscell(R.ShowProgress)
        if isnumeric(R.ShowProgress)
            R.ShowProgress=arrayfun(@num2str,R.ShowProgress,'uniformoutput',false);
        else
            R.ShowProgress=num2cell(R.ShowProgress);
        end
    end
end

%INITIALIZE PROGRESS FIGURE....
if ~isempty(R.ShowProgress)
    if islogical(R.ShowProgress)
        progressmsg='';
    else
        fh=figure;
        set(fh,'name','Templatematch progress','NumberTitle','off','renderer','opengl')
        hax=axes('pos',[0 0.01 0.5 0.95]);
        showimg(A);
        text(0.5,1,R.ShowProgress{1},'units','normalized','vert','bottom','fontname','courier','horiz','center')
        cc=zeros(2,Np);cc(:,isnan(R.Initialdu+R.Initialdv))=nan;
        hscatterA=mesh([R.pu(:) R.pu(:)]',[R.pv(:) R.pv(:)]',cc,'mesh','column','marker','.','markersize',7,'cdata',cc); %bizarrely much faster than scatter
        colormap autumn
        caxis([0 1])
        hax(2)=axes('pos',[0.5 0.01 0.5 0.95]);
        showimg(B); hold on
        hscatterB=copyobj(hscatterA,gca);
        %hscatterB=mesh(points(:,[1 1])'+duyo(1),points(:,[2 2])'+duyo(2),zeros(2,size(points,1)),'mesh','column','marker','.','markersize',5,'cdata',cc); %bizarrely much faster than scatter
        caxis([0 1])
        htext=text(1,1,'','units','normalized','vert','bottom','fontname','courier','horiz','right');
        text(0.5,1,R.ShowProgress{end},'units','normalized','vert','bottom','fontname','courier','horiz','center')
        linkaxes(hax,'xy');
        axis image, axis ij, drawnow
        hprogress=annotation(fh,'rectangle',[0 0 0 0.01],'color','none','facecolor','r');
    end
end
lastdraw=cputime;

if R.SuperSample==1
    resizefun=@(A,super)A;
else
    if license('test','image_toolbox') %(exist('imresize','file')>1)
        %use imresize if it is available. (requires image processing toolbox)
        resizefun=@(A,super)imresize(A,super);
    else
        if R.SuperSample>1
            R.SuperSample=2.^round(log2(R.SuperSample));
            resizefun=@(A,super)interp2(A,log2(super));
        else
            %         warning('image downsampling fallback for no image processing toolbox not implemented yet');
            %         resizefun=@(A,super)A;
            %         super=1;
            dwn=round(1./R.SuperSample);
            R.SuperSample=1./dwn;
            skipperfun=@(A,super)A(ceil(.5*dwn):dwn:end,floor(.5/super):dwn:end);
            %resizefun=@(A,super)skipperfun(filter2(fspecial('gaussian',[0 0]+round(2./super),1./super)),A);
            resizefun=@(A,super)skipperfun(filter2(ones(dwn)/dwn^2,A),super);
            
            %implement reshape based stacker for speed.
        end
    end
end

for ii=1:Np
    %select current point:
    p=[pu(ii) pv(ii)];
    SearchWidth=R.SearchWidth(min(numel(R.SearchWidth),ii))-1;
    SearchHeight=R.SearchHeight(min(numel(R.SearchHeight),ii))-1;
    TemplateWidth=R.TemplateWidth(min(numel(R.TemplateWidth),ii))-1;
    TemplateHeight=R.TemplateHeight(min(numel(R.TemplateHeight),ii))-1;
    Initialdu=R.Initialdu(min(numel(R.Initialdu),ii));
    Initialdv=R.Initialdv(min(numel(R.Initialdv),ii));
    
    % Actual pixel centre might differ from pu,pv because of rounding  
    % 
    Acenter=round(p) - mod([TemplateWidth TemplateHeight]/2,1);  % centre coordinate of template 
    Bcenter=round(p+[Initialdu Initialdv]) - mod([SearchWidth SearchHeight]/2,1); % centre coordinate of search region
    
    %what was actually used:
    pu(ii)=Acenter(1);
    pv(ii)=Acenter(2);
    Initialdu=Bcenter(1)-Acenter(1);
    Initialdv=Bcenter(2)-Acenter(2);


    
    try
        BB=B( Bcenter(2)+(-SearchHeight/2:SearchHeight/2)  ,Bcenter(1)+(-SearchWidth/2:SearchWidth/2),:);
        AA=A( Acenter(2)+(-TemplateHeight/2:TemplateHeight/2),Acenter(1)+(-TemplateWidth/2:TemplateWidth/2),:);
    catch
        %out of bounds... continue (and thus return a nan for that point)
        continue
    end
    
    
    
    AA=resizefun(mean(AA,3),R.SuperSample); %TODO: improve edge effects of super sampling. (need 2 extra pixels for cubic) 
    BB=resizefun(mean(BB,3),R.SuperSample);
    
    [C,uu,vv]=matchfun(AA,BB);
    [Cmax,mix]=max(C(:));
    [mix(1),mix(2)]=ind2sub(size(C),mix);
    
    meanAbsCorr(ii)=mean(abs(C(:))); %"noise" correlation level (we can accept that estimate even if we cannot find a good peak.)
    
    if ~(any(mix==1)||any(mix==size(C))) %do not accept any maxima on the edge of C
        %really simple/fast/crude sub pixel.  TODO: find bicubic interpolation max. (For now just super sample the imge for higher precision.)
        c=C(mix(1)+(-1:1),mix(2)+(-1:1));
        [uu,vv]=meshgrid(uu(mix(2)+(-1:1)),vv(mix(1)+(-1:1)));
        c=(c-mean(c(:)));c(c<0)=0; %simple and excellent performance for landsat test images... 
        c=c./sum(c(:)); 
        mix(2)=sum(uu(:).*c(:));
        mix(1)=sum(vv(:).*c(:));
%         %ALTERNATIVE 3x3 METHOD: http://www.mathworks.com/matlabcentral/fileexchange/26504-sub-sample-peak-fitting-2d
%         %by: Eric from HTWK Leipzig
%         pa = (c(2,1)+c(1,1)-2*c(1,2)+c(1,3)-2*c(3,2)-2*c(2,2)+c(2,3)+c(3,1)+c(3,3));
%         pb = (c(3,3)+c(1,1)-c(1,3)-c(3,1));
%         pc = (-c(1,1)+c(1,3)-c(2,1)+c(2,3)-c(3,1)+c(3,3));
%         %pd = (2*c(2,1)-c(1,1)+2*c(1,2)-c(1,3)+2*c(3,2)+5*c(2,2)+2*c(2,3)-c(3,1)-c(3,3));
%         pe = (-2*c(2,1)+c(1,1)+c(1,2)+c(1,3)+c(3,2)-2*c(2,2)-2*c(2,3)+c(3,1)+c(3,3));
%         pf = (-c(1,1)-c(1,2)-c(1,3)+c(3,1)+c(3,2)+c(3,3));
%         % (ys,xs) is subpixel shift of peak location relative to center (2,2)
%         mix = [(6*pb*pf-8*pe*pc) (6*pb*pc-8*pa*pf)]/(16*pe*pa-9*pb^2);
        %alternative method 2 (7x7?): http://www.mathworks.com/matlabcentral/fileexchange/46964-fastreg-zip/content//fastreg.m
        %looks vaguely similar to my own. 
        
        
        mix=mix([2 1])./R.SuperSample;
        du(ii)=mix(1)+Initialdu;
        dv(ii)=mix(2)+Initialdv;
        peakCorr(ii)=Cmax;
    end
    
    if ~isempty(R.ShowProgress)
        if islogical(R.ShowProgress)
            if (cputime-lastdraw)>.1||(cputime<lastdraw)||(ii==Np)
                backspaces=char(zeros(size(progressmsg))+8);
                progressmsg=[uint8((1:40)<((ii/Np)*40)).*'+' ''];
                progressmsg=sprintf('Templatematch [%s]  (%5.1f %+5.1f)',progressmsg,du(ii),dv(ii));
                fprintf('%s%s',backspaces,progressmsg)
                drawnow
                lastdraw=cputime;
            end
        else
            cc(:,ii)=min(max(peakCorr(ii)-meanAbsCorr(ii),0),1);
            set(htext,'string',sprintf('%+5.1f %+5.1f ',du(ii),dv(ii)))%,'units','normalized','vert','top')
            set(hprogress,'position',[0 0 ii/Np 0.01])
            set(fh,'name',sprintf('Templatematch %3.0f%%',ii*100/Np));
            if (cputime-lastdraw)>.3||(cputime<lastdraw)||(ii==Np)
                set(hscatterA,'cdata',cc);
                posB=[R.pu(:)+du(:) R.pv(:)+dv(:)];
                set(hscatterB,'cdata',cc,'xdata',posB(:,[1 1])','ydata',posB(:,[2 2])');
                zlim([-1.1 .1]) %critical as otherwise matlabs clipping plane will throw out points with z=0 in older versions of matlab. BUG.
                drawnow
                lastdraw=cputime;
            end
        end
    end
    
end

if ~isempty(R.ShowProgress)
    try
        if islogical(R.ShowProgress)
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

function [result,uu,vv]=phasecorr2(T,I)
%
%TODO-test: apply band-pass - (to remove SuperSampled high freqs, and
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
uu=-wkeep(2):wkeep(2);
vv=-wkeep(1):wkeep(1);
rows=mod(vv+c(1)-1,sz(1))+1;
cols=mod(uu+c(2)-1,sz(2))+1;
result=result(rows,cols);

function [C,uu,vv]=matNCC(T,B)
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
uu=-wkeep(2):wkeep(2);
vv=-wkeep(1):wkeep(1);


function [C,uu,vv]=myNCC(T,B)
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
uu=-wkeep(2):wkeep(2);
vv=-wkeep(1):wkeep(1);


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



