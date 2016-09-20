%% camera class - a distorted camera model.
%
% This class is an implementation of a distorted camera model.
% Note: Uses <http://docs.opencv.org/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html the same distorted camera model as opencv>.
%
%
% camera properties:
%  imgsz   - size of image in pixels [#rows, #columns]
%  f       - focal length in pixel units (two element vector [fx,fy])
%  c       - camera center in pixel coordinates (two element vector: [cx,cy])
%  k       - radial distortion coefficients. (six element vector: [k1-k6])
%  p       - tangential distortion coefficients (two element vector: [p1,p2])
%  xyz     - world coordinates of camera.
%  viewdir - [yaw,pitch,roll]. Yaw: rotation about z (0=looking east)
%                              Pitch: look up/down angle
%                              Roll: camera roll (horizon tilt).
%
% camera Dependent/derived properties:
%  (i.e. properties calculated from the camera parameters)
%  R         - camera rotation matrix calculated from camera view
%              direction (read only)
%  fullmodel - a 20-element vector containing all camera properties.
%              [camx,camy,camz,imgszy,imgszx,viewdiryaw,viewdirpitch,viewdirroll,fx,fy,cx,cy,k1-6,p1-2]
%
% camera Methods:
%  camera      - constructor
%  optimizecam - optimize the camera to mimimize misfit between
%                projected world coordinates and pixel coordinates.
%  project     - project world coordinates to image coordinates (3d->2d)
%  invproject  - project image coordinates to world coordinates (2d->3d)
%
%
% ImGRAFT - An image georectification and feature tracking toolbox for MATLAB
% Copyright (C) 2014 Aslak Grinsted (<www.glaciology.net glaciology.net>)

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

classdef camera
    
    properties
        xyz = [0 0 0]; %world coordinates of camera;
        imgsz = [100 100]; %size of image in pixels [#rows, #columns] Note: swapped xy
        viewdir = [0 0 0]; %view direction: yaw/pitch/roll. Yaw: rotation about z (0=looking east) | Pitch: look up/down angle | Roll: camera roll (horizon tilt)
        f = [5000 5000]; %focal length (two element vector).
        c = [50 50]; %camera center in pixel coordinates
        k = [0 0 0 0 0 0]; %k1-k6: radial distortion coefficients.
        p = [0 0]; %p1-p2: tangential distortion coefficients
    end
    
    properties (Dependent)
        % Camera rotation matrix calculated from viewdir. (read-only)
        R
        
        % All camera properties serialized into a 20 element vector.
        % fullmodel=[camx,camy,camz,imgszy,imgszx,viewdiryaw,viewdirpitch,viewdirroll,fx,fy,cx,cy,k1-6,p1-2]
        % Remark: Setting does not validate inputs!
        fullmodel
    end
    
    methods
        function cam = camera(varargin)
            %This is the camera constructor
            %
            %There are three ways of constructing a new camera object:
            %
            % 1. specifying all camera properties in the constructor:
            % cam = camera(xyz,imgsz[,viewdir,f,c,k,p])
            %
            % 2. Use the default camera and then modifying the camera properties
            % cam = camera()
            % cam.viewdir = [pi 0 0]; %look west.
            %
            % 3. provide a "<a href="matlab:help camera.fullmodel">fullmodel</a>" specifying all the camera
            %    properties in a single 20-element vector.
            %    (beware: no error checking).
            % cam = camera([1 1 0 1024 768 pi 0 0 1000 1000 512 384 0 0 0 0 0 0 0 0])
            %
            if nargin==0, return, end
            if nargin==1
                cam.fullmodel=varargin{1};
                return
            end
            if nargin<7,varargin{7}=[];end;
            cam.xyz = varargin{1};
            cam.imgsz = varargin{2};
            cam.viewdir=varargin{3};
            cam.f=varargin{4};
            cam.c=varargin{5};
            cam.k=varargin{6};
            cam.p=varargin{7};
            
            if length(cam.imgsz)<2, error('malformed image size.'); end
            cam.imgsz(3:end)=[];
            cam.f(end+1:2)=cam.f(end);
            if isempty(cam.c), cam.c=(cam.imgsz([2 1])+1)/2; end
            cam.k(end+1:6)=0;
            cam.p(end+1:2)=0;
        end
        
        function value = get.R(cam)
            %Camera rotation matrix calculated from viewdir
            C = cos( cam.viewdir ); S = sin( cam.viewdir );
            value = [S(3).*S(2).*C(1)-C(3).*S(1) , S(3).*S(2).*S(1) + C(3).*C(1) , S(3).*C(2); C(3).*S(2).*C(1) + S(3).*S(1), C(3).*S(2).*S(1) - S(3).*C(1) , C(3).*C(2); C(2).*C(1) , C(2).*S(1) , -S(2)];
            value(1:2,:)=-value(1:2,:);
        end
        
        function value = get.fullmodel(cam)
            value=[cam.xyz,cam.imgsz,cam.viewdir,cam.f,cam.c,cam.k,cam.p];
        end
        function cam = set.fullmodel(cam,value)
            cam.xyz=value(1:3); cam.imgsz=value(4:5); cam.viewdir=value(6:8); cam.f=value(9:10);
            cam.c=value(11:12); cam.k=value(13:18); cam.p=value(19:20);
        end
        
        function [uv,depth,inframe]=project(cam,xyz)
            % project the xyz world coordinates into image coordinates (uv)
            %
            % [uv,depth,inframe]=cam.project(xyz)
            %
            % Inputs:
            %    xyz: world coordinates
            %
            % Outputs:
            %    uv: pixel coordinates in image
            %    depth: view depth
            %    inframe: boolean vector containing whether each projected
            %    3d point is inside the frame.
            %
            if size(xyz,2)>3
                xyz=xyz';
            end
            xyz=bsxfun(@minus,xyz,cam.xyz);
            xyz=xyz*cam.R';
            xy=bsxfun(@rdivide,xyz(:,1:2),xyz(:,3));
            if any(cam.k~=0)||any(cam.p~=0) %TODO:optimize further
                r2=sum(xy.^2,2);
                r2(r2>4)=4;
                if any(cam.k(3:6)~=0)
                    a=(1+cam.k(1)*r2+cam.k(2)*r2.^2+cam.k(3)*r2.^3)./(1+cam.k(4)*r2+cam.k(5)*r2.^2+cam.k(6)*r2.^3);
                else
                    a=(1+cam.k(1)*r2+cam.k(2)*r2.^2+cam.k(3)*r2.^3);
                end
                xty=xy(:,1).*xy(:,2);
                xy=[a.*xy(:,1)+2*cam.p(1)*xty+cam.p(2)*(r2+2*xy(:,1).^2), a.*xy(:,2)+2*cam.p(1)*xty+cam.p(2)*(r2+2*xy(:,2).^2)];
            end
            uv=[cam.f(1)*xy(:,1)+cam.c(1) cam.f(2)*xy(:,2)+cam.c(2)];
            uv(xyz(:,3)<=0,:)=nan; 
            
            if nargout>1
                depth=xyz(:,3);
            end
            if nargout>2
                inframe=(depth>0)&(uv(:,1)>=1)&(uv(:,2)>=1)&(uv(:,1)<=cam.imgsz(2))&(uv(:,2)<=cam.imgsz(1)); %todo: additional constraint for negative k1 and r2>1. (See orthorectification example)
            end
        end
        
        function xyz=invproject(cam,uv,X,Y,Z,xy0)
            % Inverse projection from 2d to 3d coordinates. (pixel->world)
            %
            % There are 3 ways to call this method:
            %
            % 1. Return a set of world coordinates which when projected
            %    results in the pixel coordinates (uv).
            %
            %   xyz=cam.invproject(uv)
            %
            %
            % 2. Inverse projection constrained to DEM surface.
            %    (Project DEM to camera-view and use interpolation to find xyz.)
            %
            %   xyz=cam.invproject(uv,X,Y,Z)
            %
            %
            % 3. Minimize misfit between projected DEM point and UV. (least squares)
            %
            %   xyz=cam.invproject(uv,X,Y,Z,xy0)
            %
            %
            %
            % Inputs:
            %   uv: 2 column matrix with pixel coordinates.
            %   [X,Y,Z]: DEM. (X,Y expected to be consistent with output from meshgrid.)
            %   [xy0]: initial guess of some xy points on the DEM which are
            %          consistent with pixel coordinates in uv.
            %
            % Outputs:
            %   xyz: 3-column matrix with world coordinates consistent with
            %        pixel coordinates in uv.
            %
            %
            nanix=any(isnan(uv),2);
            anynans=any(nanix);
            if anynans
                uv(nanix,:)=[];
            end
            if nargin==2
                
                %first an exact calculation based on non-distorted model...
                depth=1000;
                xyz=[(uv(:,1)-cam.c(1))/cam.f(1) (uv(:,2)-cam.c(2))/cam.f(2)]*depth;
                xyz(:,3)=depth;
                xyz=xyz*cam.R;
                xyz=bsxfun(@plus,xyz,cam.xyz);
                %then tune that so that it is consistent with distorted model
                if any(cam.k~=0)||any(cam.p~=0)
                    E=[1 0 0;0 1 0]*cam.R; % perturb in these directions.
                    for ii=1:size(uv,1)
                        nxyz=@(m)xyz(ii,:)+m'*E;
                        %misfit=@(m)sum((cam.project(nxyz(m))-uv(ii,:)).^2);
                        %m=fminsearch(misfit,[0 0]); %THIS REQUIRES A TOOLBOX!!!
                        misfitlm=@(m)(cam.project(nxyz(m))-uv(ii,:))'.^2;
                        m=LMFnlsq(misfitlm,[0;0]);
                        xyz(ii,:)=nxyz(m);
                    end
                end
            else
                visible=voxelviewshed(X,Y,Z,cam.xyz);
                Z=Z./visible;
                xyz=nan(size(uv,1),3);
                if nargin<6
                    [uv0,~,inframe]=cam.project([X(visible(:)),Y(visible(:)),Z(visible(:))]);
                    uv0(:,3)=X(visible(:));
                    uv0(:,4)=Y(visible(:));
                    uv0(:,5)=Z(visible(:));
                    uv0=uv0(inframe,:);
                    if exist('scatteredInterpolant','file')>1 
                        Xscat=scatteredInterpolant(uv0(:,3),uv0(:,4),uv0(:,3));
                        Xscat.Points=uv0(:,1:2);
                        Yscat=Xscat; Yscat.Values=uv0(:,4);
                        Zscat=Xscat; Zscat.Values=uv0(:,5);
                    else
                        %fallback for older versions of matlab.
                        Xscat=TriScatteredInterp(uv0(:,3),uv0(:,4),uv0(:,3));  %#ok<REMFF1>
                        Xscat.X=uv0(:,1:2);
                        Yscat=Xscat; Yscat.V=uv0(:,4);
                        Zscat=Xscat; Zscat.V=uv0(:,5);
                    end
                    xy0=[Xscat(uv(:,1),uv(:,2)) Yscat(uv(:,1),uv(:,2)) Zscat(uv(:,1),uv(:,2))];
                    xyz=xy0;

                    if anynans
                        xyz(find(~nanix),:)=xyz; %find necessary because it ensures that xyz can grow.
                        xyz(find(nanix),:)=nan;
                    end
                    return
                end
                
                if Y(2,2)<Y(1,1)
                    X=flipud(X);Y=flipud(Y);Z=flipud(Z);
                end
                if X(2,2)<X(1,1)
                    X=fliplr(X);Y=fliplr(Y);Z=fliplr(Z);
                end
                
                if exist('griddedInterpolant','file')>1
                    zfun=griddedInterpolant(X',Y',Z'); %TODO: improve robustness.
                else
                    %fallback for older versions of matlab. slower
                    zfun=@(x,y)interp2(X,Y,Z,x,y);
                end
                for ii=1:length(uv)
                    %misfit=@(xy)sum((cam.project([xy zfun(xy(1),xy(2))])-uv(ii,1:2)).^2);
                    misfitlm=@(xy)(cam.project([xy(:)' zfun(xy(1),xy(2))])-uv(ii,1:2))'.^2;
                    try
                        %[xyz(ii,1:2),err]=fminunc(misfit,xy0(ii,1:2),optimset('LargeScale','off','Display','off','TolFun',0.001)); %TODO: remove dependency. can i use LMFnlsq?
                        xyz(ii,1:2)=LMFnlsq(misfitlm,xy0(ii,1:2));
                        xyz(ii,3)=zfun(xyz(ii,1),xyz(ii,2));
                        if sum(misfitlm(xyz(ii,1:2)))>2^2
                            xyz(ii,:)=nan; %do not accept greater than 2 pixel error.
                        end
                    catch
                    end
                end
                
            end
            if anynans                
                xyz(find(~nanix),:)=xyz; %find necessary because it ensures that xyz can grow.
                xyz(find(nanix),:)=nan;
            end
        end
        
        
        function [result,rmse,AIC]=optimizecam(cam,xyz,uv,freeparams)
            % Tune the camera so that projecting xyz results in uv (least squares)
            %
            % [newcamera,rmse,AIC]=cam.optimizecam(xyz,uv,freeparams)
            %
            %
            % If uv has three columns then the third column is interpreted as a
            % weight in the misfit function.
            %
            %
            % INPUTS:
            %    xyz: world coordinates.
            %    uv: target pixel coordinates.
            %        [optional 3rd column may specify weights]
            %    freeparams: a 20-element vector describing which camera
            %         parameters should be optimized. Follows same order as
            %         cam.fullmodel.
            %
            % OUTPUTS:
            %   newcamera: the optimized camera.
            %   rmse: root-mean-square-error
            %   aic: Akaike information criterion which can be used to
            %       help in determining an appropriate degree of complexity
            %       for the camera model (i.e. avoiding overfitting). [only
            %       strictly applicable for unweighted fitting]
            %
            % EXAMPLE:
            % %optimize the three view direction parameters.
            % %viewdir are the 6-8 columns in the <a href="matlab:help camera.fullmodel">camera.fullmodel</a>.
            %   [newcamera,rmse,AIC] = cam.optimizecam(xyz,uv,'00000111000000000000')
            %
            %
            
            nanrows=any(isnan(xyz),2)|any(isnan(uv),2);
            xyz(nanrows,:)=[];
            uv(nanrows,:)=[];
            
            fullmodel0=cam.fullmodel; %this describes the initial camera that is being perturbed.
            
            freeparams=~(freeparams(:)==0|freeparams(:)=='0')'; %convert to bools
            paramix=find(freeparams);
            Nfree=length(paramix);
            mbest=zeros(1,Nfree);
            
            newcam=@(m)camera( fullmodel0 + sparse(ones(1,Nfree),paramix,m,1,length(fullmodel0)) );
            
            if size(uv,2)==3
                misfit=@(m)reshape((project(newcam(m),xyz)-uv(:,1:2)).*uv(:,[3 3]),[],1);%weighted least squares
            else
                misfit=@(m)reshape(project(newcam(m),xyz)-uv,[],1);
            end
            if isnan(misfit(mbest))
                error('All GCPs must be infront of the initial camera location for optimizecam to work.'); %TODO: write better explanation. and remove requirement. 
            end
            
            [mbest,RSS]=LMFnlsq(misfit,mbest); %WORKS SUPER FAST
            
            
            Nuv=size(uv,1);
            rmse=sqrt(RSS/Nuv);
            AIC=numel(uv)*log(RSS/numel(uv)) + 2*Nfree;
            result=newcam(mbest);
        end
        
    end % methods
end % classdef