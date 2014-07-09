function targetdir=downloadDemoData(package)
%
% This code will download additional demo data from the imgraft.glaciology.net website
% provided that it has not already been downloaded.
%

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
        targetdir=fullfile(path,'cias_data')
    case {'practise' 'schneefernerkopf' 'schneeferner'}
        url='http://imgraft.glaciology.net/documentation/examples/schneefernerkopf-example/practise-schneeferner-example.zip';
        targetdir=fullfile(path,'practise_data')
    case {'imcorr' 'bindschadler' 'icestreamd'}
        url='http://imgraft.glaciology.net/documentation/examples/imcorr-examples/imcorr-exampledata.zip';
        targetdir=fullfile(path,'imcorr_data')
    otherwise
        error('unknown package')
end
       
fname=fullfile(path,'_temp.zip');
        
if length(dir(targetdir))<3
    s=urlread(url); %matlab's urlread/urlwrite does not comply with the "moved" header. Forced to do this hack... (Note this may be fragile)
    url=regexpi(s,'href="([^\n"]+)"','once','tokens');
    url=url{1};
    urlwrite(url,fname);
    unzip(fname,targetdir);
    delete(fname);
end
