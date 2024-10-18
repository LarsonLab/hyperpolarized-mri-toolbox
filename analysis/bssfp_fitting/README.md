# Pharmacokinetic Modeling and Fitting for MS-bSSFP Acquisitions

An extended, flexible pharmacokinetic model that allows for choice of acquisition type for each metabolite: 3D bSSFP, 2D GRE or 3D GRE.

Requirements:

Optimization Toolbox

Functions:

`multisite_bSSFP_fit()` - fitting function

`sim_multisite_bSSFP()` - function to simulate signals

Example scripts:

`sim_multisite_bSSFP_ex.m` - example generating simulated signals

`multisite_bSSFP_fit_ex_avgsignal.m` - example fitting kPL for an average signal in a ROI

`multisite_bSSFP_fit_ex_voxelwise.m` - example fitting kPL and lactate T2 maps

Citation: 

Sahin S, Garnæs MF, Bennett A, et al. A pharmacokinetic model for hyperpolarized 13C-pyruvate MRI when using metabolite-specific bSSFP sequences. Magn Reson Med. 2024; 92(4): 1698-1713. doi: 10.1002/mrm.30142 https://doi.org/10.1002/mrm.30142
