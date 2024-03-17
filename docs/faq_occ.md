## How do i do orientation correlation.

Orientation correlation has been added to the templatematch function. It often out performs NCC, and has good performance on Landsat 7 images with SLC off.

You can also do it manually, which used to be the solution for doing it. ImGRAFT can also do orientation correlation if you apply a simple "orientation" pre-processing step of the images.  Here's how you might implement that:

```matlab
forient=@(A)exp(i*atan2(imfilter(A,[1 0 -1],'replicate'),imfilter(A,[1;0;-1],'replicate'))); 

oA=forient(A);
oB=forient(B);

[du,dv]=templatematch(oA,oB,'method','CCF')
```

Reference: 
[Fitch et al. 2002](http://www.bmva.org/bmvc/2002/papers/95/full_95.pdf)


CCF-O