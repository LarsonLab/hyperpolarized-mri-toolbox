# Reconstruction Methods for Hyperpolarized MRI

This folder contains tools supporting reconstruction of hyperpolarized data, both MRI as well as MRS and MRSI.

## Phase Correction

`find_phase_corr.m` attempts to estimate a zero-order phase correction to put all the signal into the real channel of complex-valued data

`phase_correct_image.m` performs voxel wise phase corrections across an image

## Coil Combination

This contains code for combining multi-coil data to improve SNR, including the "RefPeak" method for estimating coil sensitivities based on the data.

[Example of Usage in the EPI demo notebooks](../demo_notebooks/EPI%20Reconstruction%20Demo.ipynb)

Reference
```
Zihan Zhu, Xucheng Zhu, Michael A. Ohliger, Shuyu Tang, Peng Cao, Lucas Carvajal, Adam W. Autry, Yan Li, John Kurhanewicz, Susan Chang, Rahul Aggarwal, Pamela Munster, Duan Xu, Peder E.Z. Larson, Daniel B. Vigneron, Jeremy W. Gordon,
Coil combination methods for multi-channel hyperpolarized 13C imaging data from human studies,
Journal of Magnetic Resonance, Volume 301, 2019, Pages 73-79, ISSN 1090-7807, https://doi.org/10.1016/j.jmr.2019.01.015.
```
## Denoising

Included here is a tensor-based denoising method for MRSI data.
```
PET by MRI: Glucose Imaging by 13C-MRS without Dynamic Nuclear Polarization by Noise Suppression through Tensor Decomposition Rank Reduction
Jeffrey R. Brender, Shun Kishimoto, Hellmut Merkle, Galen Reed, Ralph E. Hurd, Albert P. Chen, Jan Henrik Ardenkjaer-Larsen, Jeeva Munasinghe, Keita Saito, Tomohiro Seki, Nobu Oshima, Kazu Yamamoto, Peter L. Choyke, James Mitchell, Murali C. Krishna
bioRxiv 265793; doi: https://doi.org/10.1101/265793
```
For MRI data, please see the higher-order SVD Denoising method here https://github.com/UCSF-HMTRC/hp13c_EPI-hosvd_denoising
```
Kim, Y, Chen, H-Y, Autry, AW, et al. Denoising of hyperpolarized 13C MR images of the human brain using patch-based higher-order singular value decomposition. Magn Reson Med. 2021; 86: 2497–2511. https://doi.org/10.1002/mrm.28887
```


## EPI and EPSI Demos

These demos were created to accompanying the reference below, and also include Jupyter notebook versions in `/demo Notebooks/`

[EPI demo notebooks](../demo_notebooks/EPI%20Reconstruction%20Demo.ipynb)

Note that the EPSI demo including the compressed sensing does not work due to missing functions.

```
Crane, JC, Gordon, JW, Chen, H-Y, et al. Hyperpolarized 13C MRI data acquisition and analysis in prostate and brain at University of California, San Francisco. NMR in Biomedicine. 2021; 34:e4280. https://doi.org/10.1002/nbm.4280
```

## Low Rank Methods - Local Low Rank pluse Sparse and Low Rank Hankel

See Readme files in individual folders for more information

## Symmetric EPSI

Non-uniform Fourier Transform based methods to reconstruct EPSI data.

Code from reference
```
Jiang, W., Lustig, M. and Larson, P.E.Z. (2016), Concentric rings K-space trajectory for hyperpolarized 13C MR spectroscopic imaging. Magn. Reson. Med., 75: 19-31. https://doi.org/10.1002/mrm.25577
```
