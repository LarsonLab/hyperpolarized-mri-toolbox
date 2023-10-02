# Utilities
Support functions for visualizations are described below.
Further usage examples within the toolbox are linked under
the applicable function.
<p>&nbsp;</p>

## imagescn
Function to display a series of images in a single figure.
Builds upon the functionality of built-in MATLAB's
[imagesc()](https://www.mathworks.com/help/matlab/ref/imagesc.html), a function for displaying images with
scaled colors.

### Usage
[toolbox example](https://github.com/LarsonLab/hyperpolarized-mri-toolbox/blob/34536b21f67cb074ae6cfc335f32c61c278a206a/reconstruction/Local%20Low%20Rank%20plus%20Sparse/Retrospective_Simulations/t2mapping_urea_LplusS_pca.m#L523)
``` matlab
imagescn(image,scale,[rows cols]);
```

<p>&nbsp;</p>

## imagescn_overlay
Overlay a series of images onto a series of base images
representing two separate sets of information (i.e. display
hyperpolarized carbon images on proton reference images).

### Usage
``` matlab
imagescn_overlay(base_im, basescale, overlay_im, scale, [rows cols], threshold, transparency, overlay_colormap, FigureWidth, colorbarFlag);
```

<p>&nbsp;</p>

## shadedErrorBar
Creates a continuous error bar around a line plot. Supports
both symmetric and asymmetric error bar definitions.

### Usage
An example use case of this function can be found in this
repository within the 
[HP_montecarlo_evaluation](https://github.com/LarsonLab/hyperpolarized-mri-toolbox/blob/4318de2df677dbd1e6c9aed10acc97d4cf568c1e/simulations/HP_montecarlo_evaluation.m#L446) framework.