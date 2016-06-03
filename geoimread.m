function [A,x,y,I]=geoimread(filename,varargin)
%GEOIMREAD reads a sub region of a geotiff or geojp2 image.
%
%
%% Syntax
%
% A = geoimread(filename)
% A = geoimread(filename,xlim,ylim)
% A = geoimread(filename,latlim,lonlim)
% A = geoimread(...,buffer)
% [A,x,y,I] = geoimread(...)
% geoimread(...)
%
%
%% Description
%
% A = geoimread(filename) returns the full image given by a filename. This
% syntax is equivalent to A = geotiffread(filename).
%
% A = geoimread(filename,xlim,ylim) limits the region of the geotiff file to
% the limits given by xlim and ylim, which are map units (usually meters) relative
% to the data projection. For example, if the geotiff is projected in Texas Centric
% Mapping System/Lambert Conformal coordinates, xlim and ylim will have units of
% meters relative to the origin (100 W, 18 N). xlim and ylim can be multimensional,
% in which case the limits of the map will be taken as the limits of the limits of
% the distribution of all points in xlim, ylim.
%
% A = geoimread(filename,latlim,lonlim) if no values in xlim, ylim exceed
% normal values of latitudes and longitudes, geoimread assumes you've entered
% limits in geographic coordinates of degrees latitude and longitude. The first
% input is latitude, the second input is longitude.
%
% A = geoimread(...,buffer) adds a buffer in map units (usually meters or feet) to the
% limits of the region of interest.  This may be useful if you want to load an image
% surrounding scattered lat/lon data.  If you'd like an extra 2 kilometers of image around
% your data, enter 2000 as the buffer.  If buffer is a two-element vector, the first
% element is applied to the left and right extents of the image, and the second element
% is applied to the top and bottom extents of the image.
%
% [A,x,y,I] = geoimread(...) also returns pixel center coordinates (x,y) of the
% output image and a geotiff info structure I. I is a useful input for projfwd and projinv.
%
% geoimread(...) without any outputs shows the output image A without loading
% any data into the workspace.
%
%
%% Examples:
%
% % Show a whole geotiff:
% geoimread('boston.tif');
%
% % Compare results from above to a subset geotiff:
% mapx = [765884 766035 766963]; % units are feet
% mapy = [2959218 2957723 2958972];
% geoimread('boston.tif',mapx,mapy)
%
% % Or if you have coordinates in lat/lon and you want a 500 foot buffer:
% lat = [42.3675288 42.3634246 42.3668397];
% lon = [-71.0940009 -71.0934685 -71.0900125];
%
% geoimread('boston.tif',lat,lon,500);
%
%% Author Info:
%
% (c) Aslak Grinsted 2014- (http://www.glaciology.net/)
%     & Chad A. Greene (http://chadagreene.com/)
%
%%
%
% See also GEOTIFFREAD, GEOTIFFINFO, PIXCENTERS, and PROJFWD.

%%  CHAD'S CHANGES:
% The following changes were made by Chad A. Greene (http://chadagreene.com/)
% of the University of Texas Institute for Geophysics (UTIG) on Sept 17, 2014:
%
% * More input checking and error messages.
%
% * Clarified syntax and description in header.
%
% * The fileparts line now only writes the file extension because other outputs went unused.
%
% * If geographic coordinate limits are inferred, they are now ordered latlim,lonlim. <-- **FUNCTIONALITY CHANGE**
%
% * Limits xlim, ylim or latlim, lonlim can be scalar, vector, or matrix-- xlim is now taken
%   as xlim = [min(xlim(:)) max(xlim(:))]. This will save a small step if you have some data points
%   given by x,y.  Now you can simply enter your data coordinates and geoimread will
%   figure out the limits.
%
% * Output variable I now has correct corner coordinates for the subset image instead of
%   the previous version which returned corner coordinates of the full original image.
%
% * A buffer can be specified around input limits.
%
%
% Syntax before the changes:
% geoimread('myimage.tif',[min(easting)-buffer_m max(easting)+buffer_m],...
%     [min(northing)-buffer_m max(northing)+buffer_m]);
%
% Syntax after the changes:
% geoimread('myimage.tif',easting,northing,buffer_m)
%
%

%TODO: support downsampling (ReductionLevel parameter in imread)
%TODO: use map2pix and latlon2pix instead of projfwd and pixcenters. more robust if it is a rotational coordinate system.


%% Set defaults:

usegeocoordinates = false;
returnNONsubsetImage = true;
buffer_x = 0;
buffer_y = 0;

%% Input checks:

% Check for mapping toolbox:
% assert(license('test','map_toolbox')==1,'geoimread requires Matlab''s Mapping Toolbox.')

% Check file type:
assert(isnumeric(filename)==0,'Input filename must be a string.')
[~,~,ext] = fileparts(filename);
switch upper(ext)
    case {'.JP2' '.JPEG2000' '.GEOJP2'}
        I = jp2tiffinfo(filename);
    case {'.TIF' '.TIFF' '.GTIF' '.GTIFF'}
        I = robustgeotiffinfo(filename);
    otherwise
        error('Unrecognized image file type. Must be tif, tiff, gtif, gtiff, jp2, jpeg2000, or geojp2.')
end


% Parse optional inputs:
if nargin>1
    returnNONsubsetImage = false;
    assert(nargin>2,'If you specify an xlim or latlim, you must specify a corresponding ylim or lonlim.')
    
    % Parse limits:
    xlimOrLatlim = varargin{1}(:);
    ylimOrLonlim = varargin{2}(:);
    
    assert(isnumeric(xlimOrLatlim)==1,'Input xlim or latlim must be numeric.')
    assert(isnumeric(ylimOrLonlim)==1,'Input ylim or lonlim must be numeric.')
    
    % Assume geo coordinates if no input limits exceed normal lat/lon values:
    if max(abs(xlimOrLatlim))<=90 && max(abs(ylimOrLonlim))<=360
        usegeocoordinates = true;
    end
    
    % Parse buffer:
    if nargin>3
        buffer_m = varargin{3};
        assert(isnumeric(buffer_m)==1,'Buffer value must be either a scalar or a two-element array.')
        assert(numel(buffer_m)<3,'Buffer value must be either a scalar or a two-element array.')
        buffer_x = buffer_m(1);
        if isscalar(buffer_m)
            buffer_y = buffer_m(1);
        else
            buffer_y = buffer_m(2);
        end
    end
    
    if nargin>4
        error('Too many inputs in geoimread.')
    end
    
end


%% Begin work:

% Get pixel coordinates of full (non-subset) image:
[x,y]=robustpixcenters(I);

if nargout==0
    fprintf('No outputs specified.\n')
    fprintf('X-range: %.1f:%.1f:%.1f\n',x(1),x(2)-x(1),x(end))
    fprintf('Y-range: %.1f:%.1f:%.1f\n',y(1),y(2)-y(1),y(end))
    clear A x y I
    return
end

% Set xlim and ylim depending on user inputs:
if returnNONsubsetImage
    xlimOrLatlim = x(:);
    ylimOrLonlim = y(:); 
end

if usegeocoordinates
    % lat/lon limits switch to x/y limits here:
    if ~strcmp(I.ModelType,'ModelTypeGeographic')
        assert(license('test','map_toolbox')==1,'Mapping toolbox needed to project between lat/lon limits and x,y limits. Specify limits in x,y coordinates.')
        [xlimOrLatlim,ylimOrLonlim]=projfwd(I,xlimOrLatlim,ylimOrLonlim);
    end
end


xlim = [min(xlimOrLatlim)-buffer_x max(xlimOrLatlim)+buffer_x];
ylim = [min(ylimOrLonlim)-buffer_y max(ylimOrLonlim)+buffer_y];

 
% Rows and columns of pixels to read:
rows=find((y>=ylim(1))&(y<=ylim(2)));
cols=find((x>=xlim(1))&(x<=xlim(2)));




%% Display messages if region of interest is partly or wholly outside the image region:

if xlim(1)<min(x)||xlim(2)>max(x)
    disp('geoimread limits extend beyond the available image output in the x direction.')
end

if ylim(1)<min(y)||ylim(2)>max(y)
    disp('geoimread limits extend beyond the available image output in the y direction.')
end

if isempty(rows)||isempty(cols)
    error('No image coordinates can be found inside the specified limits.')
end

%% Load region of interest:
reductionlevel=0;
if reductionlevel==0
    rows=sort(rows([1 end]));
    cols=sort(cols([1 end]));
else
    %% Load region of interest:
    dpx=2^reductionlevel;
    rows=round(rows/dpx);cols=round(cols/dpx);
    
    rows=sort(rows([1 end]));
    cols=sort(cols([1 end]));
end
x=x(cols(1):cols(end));
y=y(rows(1):rows(end));

A=imread(filename,'PixelRegion',{rows cols});


%% Update info structure to more accurately reflect the new image:

if nargout == 4
    I.FileSize = numel(A);
    I.Height = size(A,1);
    I.Width = size(A,2);
    try
        I.TiePoints.WorldPoints.X = x(1);
        I.TiePoints.WorldPoints.Y = y(1);
        I.SpatialRef.RasterSize = [size(A,1),size(A,2)];
        I.RefMatrix(3,1) = x(1);
        I.RefMatrix(3,2) = y(1);
        I.BoundingBox = [min(x) min(y); max(x) max(y)];
        I.CornerCoords.X = [min(x) max(x) max(x) min(x)];
        I.CornerCoords.Y = [max(y) max(y) min(y) min(y)];
        %TODO: check whether GTRasterTypeGeoKey is RasterPixelIsArea or RasterPixelIsPoint
        I.CornerCoords.Row = .5 + [0 0 size(A,1) size(A,1)]; %TODO: is this .5 always true?
        I.CornerCoords.Col = .5 + [0 size(A,2) size(A,2) 0];
        [I.CornerCoords.Lat,I.CornerCoords.Lon] = projinv(I,I.CornerCoords.X,I.CornerCoords.Y);
        I.GeoTIFFTags.ModelTiepointTag(4) = x(1);
        I.GeoTIFFTags.ModelTiepointTag(5) = y(1);
        I.SpatialRef.XLimWorld = [min(x),max(x)];
        I.SpatialRef.YLimWorld = [min(y),max(y)];
    catch,end
end

%% Clean up:






function I=jp2tiffinfo(fname)
%
% This is a function that attempts to extract the geo-info from a jp2 file.
% There are two "standards" for how this can be stored:
% * GeoJP2: an embedded geotiff is stored in a "uuid" box with a specific id.
% * GMLJP2: an embedded gml file is in an asoc box.
% (extremely annoying that there are multiple "standards".)
%

%Documents:
% JP2 standard: http://www.jpeg.org/public/15444-1annexi.pdf
% GEOJP2 standard: http://wiki.opf-labs.org/download/attachments/11337762/15444-1annexi.pdf
% GMLJP2 standard: http://www.opengeospatial.org/standards/gmljp2

fid=fopen(fname,'r','ieee-be');
lblcontents='';
while ~feof(fid)
    lbox=fread(fid,1,'uint32');
    type=fread(fid,[1 4],'uint8=>char');
    switch lbox
        case 1
            lbox=fread(fid,1,'uint64');lbox=lbox-16;
        case 2:7 %"reserved reserved for ISO use" (????)
            lbox=8;%No idea why ... there is no info on how to interpret values 2-7. 
            %Setting this to 8 works for sentinel-2 images. 
        otherwise
            lbox=lbox-8;
    end
    
    switch type
        case 'asoc'
            %the asoc box is a container with other boxes. Set lbox to zero
            %so that we also parse the containing boxes. 
            lbox=0;
        case 'lbl '
            %in gmljp2 the xml box is preceded by a lbl box. Read the
            %lbl box contents for later...
            lblcontents=fread(fid,[1 lbox],'uint8=>char');
            lbox=0;
        case 'xml '
            if (strfind(lblcontents,'gml')>0) %gmljp2:
                gmlcontents=fread(fid,[1 lbox],'uint8=>char');
                fout=fopen('.temp.xml','w');
                fwrite(fout,gmlcontents);
                fclose(fout);
                X=xmlread('.temp.xml');
                delete('.temp.xml');
                
                Elements=X.getElementsByTagName('gml:origin');
                origin=str2num(Elements.item(0).getTextContent());
                Elements=X.getElementsByTagName('gml:offsetVector');
                R=[];
                for ii=1:Elements.getLength
                    R{ii}=str2num(Elements.item(ii-1).getTextContent());
                end
                R=cell2mat(R');
                %ModelPixelScaleTag: [15 15 0]
                %ModelTiepointTag: [0 0 0 4.695e+05 7240500 0]
                assert((R(1,2)==0)&(R(2,1)==0),'unexpected offsetvector in GML');
                I=imfinfo(fname);
                I.GML=gmlcontents;
                I.ModelPixelScaleTag=[R(1,1) -R(2,2) 0]; %order and sign..
                I.ModelTiepointTag=[0 0 0 origin 0];
                fclose(fid);
                return %we have what we need - return
            end
        case 'uuid' %
            uuid=fread(fid,[1 16],'uint8=>char');
            lbox=lbox-16;
            geo=[177 75 248 189 8 61 75 67 165 174 140 215 213 166 206 3];
            if all(uuid==geo)
                fout=fopen('.temp.tif','w');
                contents=fread(fid,lbox,'uint8');lbox=0;
                fwrite(fout,contents);
                fclose(fout);
                fclose(fid);
                I=geotiffinfo('.temp.tif');
                m=imfinfo(fname); % a little silly to use imfinfo when I already have a tag reader
                I.Height=m.Height;
                I.Width=m.Width;
                delete('.temp.tif');
                
                fclose(fid);
                return %we have what we need - return
            end
    end

    fseek(fid,lbox,0);
    
end
fclose(fid);



%--- BELOW is to make geoimread work even if no mapping toolbox (but with reduced functionality) ---

function I = robustgeotiffinfo(fname)
if license('test','map_toolbox')
    I=geotiffinfo(fname);
else
    I=imfinfo(fname);
    % I.ModelType='ModelTypeProjected';  %TODO: fix
    %     %TODO: generate home-made refmatrix(?)....
    %     if isfield(tags, 'ModelTransformationTag') && numel(tags.ModelTransformationTag) >= 8
    %         geoimread does not work for rotated systems
    %
    %     else %use ModelPixelScaleTag instead
    %         dx =  I.ModelPixelScaleTag(1); dy = -I.ModelPixelScaleTag(2);
    %         x0 = I.ModelTiepointTag(4) - dx * I.ModelTiepointTag(1);
    %         y0 = I.ModelTiepointTag(5) - dy * I.ModelTiepointTag(2);
    %         J = [dx 0; 0 dy];
    %     end
    %     I.RefMatrix=[flipud(J); x0-J(1,1)-J(1,2)];
end

function [x,y]=robustpixcenters(I)
if license('test','map_toolbox') && isfield(I,'GeoTIFFCodes')
    [x,y]=pixcenters(I);
else
    %I have not read documentation... but this only works for rectilinear systems.
    assert(I.ModelPixelScaleTag(3)==0,'unexpected ModelPixelScaleTag format.');
    assert(all(I.ModelTiepointTag(1:3)==0),'unexpected ModelTiepointTag format.');
    x=((0:I.Width-1)-I.ModelTiepointTag(1))*I.ModelPixelScaleTag(1)+I.ModelTiepointTag(4);
    y=((0:I.Height-1)-I.ModelTiepointTag(2))*-I.ModelPixelScaleTag(2)+I.ModelTiepointTag(5);
end

