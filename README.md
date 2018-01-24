The Fast Iterative Digital Volume Correlation Algorithm (FIDVC) is the next generation DVC algorithm providing significantly improved signal-to-noise, and large (finite) deformation (incl. large rotations and image stretches) capture capability at low computational cost (please see [Bar-Kochba, Toyjanova et al., Exp. Mechanics, 2014](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst) for more details).

### Important pages
* [Download latest version v1.2.4!](https://github.com/FranckLab/FIDVC/releases)
* [Example data](https://drive.google.com/folderview?id=0ByhZqlrbo5srSmU2ZW1TOXpfVkE&usp=sharing)
* [FAQ](https://github.com/FranckLab/FIDVC#faq)
* [Questions/Issues](https://github.com/FranckLab/FIDVC/issues)
* [Bug Fixes/history](https://github.com/FranckLab/FIDVC/wiki/Bug-Fixes!)
* [Franck Lab](http://franck.engin.brown.edu)
 
## Purpose
The following implementation contains the Matlab m-files for our FIDVC algorithm. The FIDVC algorithm determines the 3D displacement fields between consecutive volumetric image stacks. To provide maximum user versatility the current version is the CPU-based version, which runs a little bit slower (~ 5 - 15 minutes/per stack) than the GPU-based implementation. As Matlab’s GPU-based subroutines become more efficient we hope to provide the GPU-based version at a later release date. 

## Running FIDVC

### C Compiler
To run you need a compatible C compiler. Please see
(http://www.mathworks.com/support/compilers/R2015a/index.html)

### Input 3D Image Stack Requirements
* To check if the 3D image stack have the required speckle pattern and intensity values for correlation please use our [DVC simulator](https://github.com/FranckLab/DVC-Simulator).
* The 3D image stack need to be saved in a 3 dimensional matrix (intensity values are stored at x, y and z position) in **vol*.mat** files.  
* We recommend that the input image stack at each dimension should have at least 1.5 times of the subset size as the number of pixels. The default subset size is 128x128x64, so we recommend that the minimum input volume size should be 192x192x96.
* The size of the input image stack should be divisible by 0.5 times the size of the subset. 

### Running included example case
1. Make sure that the main files and the supplemental m files (from file exchange) are added to the path in Matlab.
2. Download and save the [example volume data](https://drive.google.com/folderview?id=0ByhZqlrbo5srSmU2ZW1TOXpfVkE&usp=sharing) in the example folder. 
3. Run the exampleRunFile.m file to and compare its displacement outputs to the contour plots in the referenced paper ([Bar-Kochba, Toyjanova et al., Exp. Mechanics, 2014](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst))

### Health warning!
FIDVC requires a 3D stack to be read in, which depending on the volume size can require a **large amount of RAM** in Matlab.

## Files
* Main files
 - addDisplacements.m
 - checkConvergenceSSD.m
 - DVC.m
 - filterDisplacements.m
 - funIDVC.m
 - IDVC.m
 - removeOutliers.m
 - volumeMapping.m

* Supplement m files from the MATLAB file exchange:
 - gridfit.m
 - inpaint_nans.m
 - inpaint_nans3.m
 - mirt3D_mexinterp.cpp
 - mirt3D_mexinterp.m
 - mirt3D_mexinterp.mexw64

* Example Run files
 - exampleRunFile.m
 - [example volume data](https://drive.google.com/folderview?id=0ByhZqlrbo5srSmU2ZW1TOXpfVkE&usp=sharing) (vol00.mat, vol01.mat, resultsFIDVC.mat, outputREsults.pdf, matlab_workspace.mat).

## FAQ
**What are the requirements for the input 3D image stack?**

Please refer to [input 3D Image Stack Requirements](https://github.com/FranckLab/FIDVC#input-3d-image-stack-requirements).

**Can I use FIDVC for finding displacement fields in 2D images?**

No. But you can use [FIDIC](https://github.com/FranckLab/FIDIC), this is 2D version of FIDVC for finding displacments in 2D images. 

## Cite
If used please cite:
[Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast iterative digital volume correlation algorithm for large deformations. Experimental Mechanics. doi: 10.1007/s11340-014-9874-2](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst)

```bibtex
@article{bar2014fast,
  title={A fast iterative digital volume correlation algorithm for large deformations},
  author={Bar-Kochba, E and Toyjanova, J and Andrews, E and Kim, K-S and Franck, C},
  journal={Experimental Mechanics},
  pages={1--14},
  year={2014},
  publisher={Springer}
}
```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/FIDVC#faq) and [Questions/Issues](https://github.com/FranckLab/FIDVC/issues). Add a new question if similar issue hasn't been reported. We shall help you at the earliest. The author's contact information can be found at [Franck Lab](http://franck.engin.brown.edu).
