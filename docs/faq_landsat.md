## How do i get velocities from Landsat images?

Here are a few tips:

You can download landsat scenes from Earth Explorer. Usually it is best to use the high-resolution pan-chromatic band in feature tracking. 
It is important that both images are taken from the same view point as the georectification in the L1G product is not 100% exact. However, if you choose an image pair where both images share the same viewpoint then these georectification errors will tend to cancel out. You should therefore only track features in images that share the same row and path number.

Landsat scenes are rather large, and you may run into memory issues if you attempt to track whole scenes at once. I have written a small function that lets you load a smaller study area of a large geotiff. It is called geoimread and can be downloaded on matlabcentral. 

It is generally a good idea to high-pass filter the input images as a pre-processing step. See this FAQ on preprocessing for how to do that.

Once you have loaded the images, and applied any pre-filtering, then you can just call templatematch to get the displacement. Here you can use these two examples as templates: Batura, Bindschadler.

A 'default' template size that usually works is around 21 pixels on a side. But larger templates may work better in slow moving regions. Larger templates have the disadvantage that they are more sensitive to shear. So, large templates work best when ice flow is relatively uniform on the spatial scales of the template. 

In areas with slow ice flow, you need larger temporal baselines in the image pairs to accurately estimate movement. 

It can be a good idea to make an animated gif of the two images in the image pair. It is a great way to get an intuition about what is going on. It is great for trouble shooting. 

