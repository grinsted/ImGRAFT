## Tip: Batch processing of a time lapse sequence

When batch processing a time-lapse sequence then it is a good idea to separate the work flow into smaller digestible chunks and then save the results along the way. This outline shows how you might construct such a script.  

* Determine camera orientation and internal parameters of a master image.
* Get as many GCPs for that image as possible.
* Use camera.optimizecam to determine unknown parameters. (see Schneefernerkopf for an example).
* Determine the camera orientation of all the remaining images with respect to master image. This can be done by tracking the features of background movement. See the [Engabreen](demoengabreen.md) demo for an example.
* Pick a set of real world coordinates to feature track. See the [Engabreen](demoengabreen.md) demo for an example.
* Determine the pixel displacement of all candidate image pairs.
  - Some proposed criteria that can be used to select candidate pairs:
  - Roughly same time of day (to minimize issues from shadow movement). E.g. dt>1day and dt<14days.
  - Criteria that removes ill-suited images: under-exposed, cloudy, ice/water on lens. 
  - This can be done using manually or automatically using e.g. image statistics, or exposure settings. 
* For each pair: Convert all pixel locations of tracked features to 3d positions, and then to velocities. 

See [Engabreen](demoengabreen.md) for an example.

I would script each of these steps separately and save the results along the way. 

    
Many of these steps require you to do some task for each file. Here's how you would typically solve that in Matlab:

```matlab
files=dir(fullfile(imagefolder,'img_*.jpg'));
for ii=1:length(files)
 fname=fullfile(imagefolder,files(ii).name);

 [... insert your code here ...]

end
```

## Refinements:
There can be a slight error in calculating the change in view direction of the image with respect to the master image. Especially if the light and/or snow conditions are very different between images. This introduces a small error in the 2d->3d coordinates projection. You can reduce this error by having multiple master images depending on the conditions (early/late season, direct/diffuse light, morning/evening). When you calculate velocities then you are differencing two the locations derived from an image pair. If part of the location error is shared then this will cancel out. It can therefore be advantageous to derive the camera parameters for image B from camera A rather than from the master camera.

