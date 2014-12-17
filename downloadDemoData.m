function targetdir=downloadDemoData(package)
%% download additional imgraft demo data
%
% targetdir=downloadDemoData(package)
% 
% This code will download additional demo data from the imgraft.glaciology.net website
% provided that it has not already been downloaded.
%
%
% inputs:
%    Package can be: 'cias', 'practise', or 'imcorr'
%
% outputs:
%    targetdir is the location of the downloaded files.
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

if nargin==0
    downloadDemoData('cias');
    downloadDemoData('practise');
    downloadDemoData('imcorr');
    return
end

path=fileparts(mfilename('fullpath'));
path=fullfile(path,'demos');

switch lower(package)
    case {'cias' 'muragl' 'batura'}
        url='http://imgraft.glaciology.net/documentation/examples/cias-example/cias_data.zip';
        targetdir=fullfile(path,'cias_data');
        testfile=fullfile(targetdir,'batura_2000.tif');
    case {'practise' 'schneefernerkopf' 'schneeferner'}
        url='http://imgraft.glaciology.net/documentation/examples/schneefernerkopf-example/practise-schneeferner-example.zip';
        targetdir=fullfile(path,'practise_data');
        testfile=fullfile(targetdir,'dem30m.mat');
    case {'imcorr' 'bindschadler' 'icestreamd'}
        url='http://imgraft.glaciology.net/documentation/examples/imcorr-examples/imcorr-exampledata.zip';
        targetdir=fullfile(path,'imcorr_data');
        testfile=fullfile(targetdir,'conv_89.png');
    otherwise
        error('unknown download package')
end
       
fname=fullfile(path,'_temp.zip');
        
if ~exist(testfile,'file')
    s=urlread(url); %matlab's urlread/urlwrite does not comply with the "moved" header. Forced to do this hack... (Note this may be fragile)
    url=regexpi(s,'href="([^\n"]+)"','once','tokens');
    url=url{1};
    urlwrite(url,fname);
    unzip(fname,targetdir);
    delete(fname);
end
