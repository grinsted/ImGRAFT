# ImGRAFT
## An image georectification and feature tracking toolbox for MATLAB

On the [ImGRAFT](http://imgraft.glaciology.net) web site you will find more detailed information on the project such as:
* [Documentation](http://imgraft.glaciology.net/documentation)
* [Examples](http://imgraft.glaciology.net/documentation/examples)
* [Project news](http://imgraft.glaciology.net/news)

This software is open source (See licensing details elsewhere). In addition to the formal licensing terms, We would greatly appreciate an acknowledgement. Preferably in the form of a citation (once the paper is out) and a link to the web-page. 


Note: Example data not released yet because licensing terms has not yet been decided. 



## File descriptions:

* example.m
	* A working example feature tracking two images of Engabreen and calculating real world coordinates. This is a good place to start
* camera.m
	* A model of a camera with lens distortion. Used to project between 2d and 3d space. Can work with DEM data.
* templatematch.m
	* finds features in image A in image B using template matching. 
* voxelviewshed.m
	* A fast viewshed calculation of a DEM. 
* LMFnlsq.m
	* helper function used in the least squares optimization of the camera parameters. (By M.Balda)
	

	
## Licensing

The majority of the code is licensed under a very permissive MIT license, but some routines and example data is licensed under other terms.	See licensing details in LICENSE.txt and individual files. 

This software package includes the following open source codes licensed under other terms:

* LMFnlsq.m 
	* Copyright Miroslav Balda. This is an implementation of the Levenberg-Marquardt algorithm as modified by Fletcher. It is used in the least squares optimization of the camera parameters. See licensing details in LMFnlsq.m 


## [Acknowledgements](http://imgraft.glaciology.net/acknowledgements)

This software has been developed at [Centre for Ice and Climate](http://www.iceandclimate.nbi.ku.dk), Niels Bohr Institute, University of Copenhagen as part of the [SVALI project](http://www.ncoe-svali.org/). SVALI is a part of the [Top-level Research Initiative](http://www.norden.org/) (TRI), which is a major Nordic collaborative venture for studies of climate, energy and the environment. We are also grateful to Miriam Jackson and [NVE](http://nve.no) who has helped facilitate the Engabreen fieldwork and contributed with data.
 

