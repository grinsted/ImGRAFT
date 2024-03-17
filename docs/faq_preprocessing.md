## Preprocessing images for feature tracking

In some cases you may improve the feature tracking by pre-processing the images. The aim of the preprocessing may be to reduce effects from different lighting conditions, remove noise, or reduce the weight of outliers. Here are a set of proposed pre-processing steps that may improve your results. 

Many of the filters below are written using matlabs @-syntax for anonymous functions. You would apply the filters like this:

```matlab
A=im2single(imread('test.tif'));
hpass=@(A,sigma)A-imfilter(A,fspecial('gaussian',sigma*3,sigma)); %define a filter function
A=hpass(A,3); %apply the filter
imagesc(A)
```


## High-pass filtering
Feature tracking using NCC will attempt to align bright and dark regions of the image. These may correspond to illuminated vs shadow regions of the ice and thus depend strongly on the direction of the illumination. For that reason you may want to apply a high pass filter to focus on finer details than large patches of bright and dark. This can be written in matlab like this:

```matlab
hpass=@(A,sigma)A-imfilter(A,fspecial('gaussian',sigma*3,sigma))
```

## Local histogram equalization
Histogram equalization will bring in outliers so that they are not allowed to dominate the output.

```matlab
histequalfilt=@(A,tilesize)adapthisteq(A,'NumTiles',ceil(size(A)./tilesize));
```

See help on adapthisteq from the image processing toolbox. 


## Noise reduction
See e.g. the wiener2 function in the image processing toolbox. 

Combination of filters
You can combine multiple filters.  E.g.

```matlab
myfilter=@(A)histequalfilt(hpass(A,3),50);
```

## Structure texture separation
Structure texture separation separates the image into a structure component which can be thought of as the background color, and a residual called the texture. In this respect the structure component is similar to a low-pass filtered image except that structure-texture methods allow for sharp transitions between regions of different brightness. 
I have not tried structure-texture separation techniques yet but have read papers that indicate that this may be a very fruitful venue for tracking between scenes where the sun position has changed.  It could help suppress the sharp brightness change at the edge of a shadow from a mountain. 

I recommend this matlab package: http://www.cse.cuhk.edu.hk/leojia/projects/texturesep/index.html

## Gradient orientation transformation
See also Orientation correlation -- This is now a built-in method of templatematch.