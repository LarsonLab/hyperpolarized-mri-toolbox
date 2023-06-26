# RF Pulses for Hyperpolarized MRI

## Spectral-Spatial RF Pulses

Included here is a self-contained Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package.  See the `spectral_spatial/` folder for more detailed information and references.

## Variable, Metabolite-Specific Flip Angle methods

Included here are scripts to general variable, metabolite-specific flip angle schemes tailored to hyperpolarized MRI.  See the `variable_flip_angles/` folder for more detailed information and references.

## Slice Profile

In hyperpolarized MRI, uneven depletion of the hyperpolarized magnetization leads to distortions of the RF pulse slice profile after repeated RF pulses.
The code here simulates this effect and computes potential compensation factors.

### Relevant citations

Martin H. Deppe, Kevin Teh, Juan Parra-Robles, Kuan J. Lee, Jim M. Wild,
Slice profile effects in 2D slice-selective MRI of hyperpolarized nuclei,
Journal of Magnetic Resonance,
Volume 202, Issue 2,
2010,
Pages 180-189,
ISSN 1090-7807,
https://doi.org/10.1016/j.jmr.2009.11.003.

Gordon, J.W., Milshteyn, E., Marco-Rius, I., Ohliger, M., Vigneron, D.B. and Larson, P.E.Z. (2017), Mis-estimation and bias of hyperpolarized apparent diffusion coefficient measurements due to slice profile effects. Magn. Reson. Med., 78: 1087-1092. https://doi.org/10.1002/mrm.26482

Walker, CM, Gordon, JW, Xu, Z, et al. Slice profile effects on quantitative analysis of hyperpolarized pyruvate. NMR in Biomedicine. 2020; 33e4373. https://doi.org/10.1002/nbm.4373
