The Fast Iterative Digital Volume Correlation Algorithm (FIDVC) is the next generation DVC algorithm providing significantly improved signal-to-noise, and large (finite) deformation (incl. large rotations and image stretches) capture capability at low computational cost (please see [Bar-Kochba, Toyjanova et al., Exp. Mechanics, 2014](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst) for more details).

## Important pages
* [Download latest version v1.1!](https://github.com/FranckLab/FIDVC/releases)
* [Example data](http://franck.engin.brown.edu/~christianfranck/FranckLab/downloads/FIDVC.zip)
* [FAQ](https://github.com/FranckLab/FIDVC/blob/master/README.md#faq)
* [Questions/Issues](https://github.com/FranckLab/FIDVC/issues)
* [Franck Lab](http://franck.engin.brown.edu)
 
# Purpose
The following implementation contains the Matlab m-files for our FIDVC algorithm. The FIDVC algorithm determines the 3D displacement fields between consecutive volumetric image stacks. To provide maximum user versatility the current version is the CPU-based version, which runs a little bit slower (~ 5 - 15 minutes/per stack) than the GPU-based implementation. As Matlab’s GPU-based subroutines become more efficient we hope to provide the GPU-based version at a later release date. 

# Running FIDVC

### C Compiler
To run you need a compatible C compiler. Please see
(http://www.mathworks.com/support/compilers/R2015a/index.html)

## Input 3D Image Stack Requirements
* To check if the 3D image stack have the required speckle pattern and intensity values for correlation please use our [DVC simulator](https://github.com/FranckLab/DVC-Simulator).
* The 3D image stack need to be saved in a 3 dimensional matrix (intensity values are stored at x, y and z position) in **vol*.mat** files.  
* We recommend that the input image stack at each dimension should have at least 1.5 times of the subset size as the number of pixels. The default subset size is 128x128x64, so we recommend that the minimum input volume size should be 192x192x96.
* The size of the input image stack should be divisible by 0.5 times the size of the subset. 

## Steps


Main files:
addDisplacements.m
checkConvergenceSSD.m
DVC.m
filterDisplacements.m
funIDVC.m
IDVC.m
removeOutliers.m
volumeMapping.m

## Health warning!
The DVC simulator requires a 3D stack to be read in, which depending on the volume size can require a **large amount of RAM** in Matlab, and thus when running the simulator from a laptop, we strongly suggest to use the resize option within the GUI to decrease the total stack dimensions. My macbook air can only handle 128 x 128 x 96 sufficiently. 
