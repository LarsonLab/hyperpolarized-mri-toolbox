# Variable Flip Angle RF Pulse Design scripts

Original Authors:  Peder E. Z. Larson, Ann Xing, Sonam Machingal, John Maidens

These Matlab functions are for designing variable flip angle schemes for hyperpolarized MR applications.  All functions are documented within.  

## Description of functions:

### vfa_const_amp()

RF-compensated variable flip angle scheme for constant signal.
Can also include compensation for T1 or an effective T1, aiming for constant signal amplitude.  Demo in `vfa_example_script.m`

References: 
```
Zhao, L., Mulkern, R., Tseng, C.-H., Williamson, D., Patz, S., Kraft, R., Walsworth, R. L., Jolesz, F. A., and Albert, M. S. Gradient-echo imaging considerations for hyperpolarized 129xe mr. J Magn Reson B 113, 2 (November 1996), 179–183.
```
```	
Xing Y, Reed GD, Pauly JM, Kerr AB, and Larson PEZ. Optimal Variable Flip Angle Schemes For Dynamic Acquisition Of Exchanging Hyperpolarized Substrates.  J Magn Reson. 2013 Sep; 234:75-81.  http://www.ncbi.nlm.nih.gov/pubmed/23845910
```
### vfa_opt_signal()

flip angle schemes to maximize the total signal when summed across all acquisitions. 
Theoretically derived to include compensation for T1 decay.   Demo in `vfa_example_script.m`. 
Emperically derived to be optimal for metabolic products as well in `optimal_SNR_demo.m`

References:
```
Nagashima, K. Optimum pulse flip angles for multi-scan acquisition of hyperpolarized nmr and mri. J Magn Reson 190, 2 (February 2008), 183–188.
```
### mvfa_const_amp()

Constant amplitude signal flip angles in the presence of both T1 decay and metabolic conversion.  Demo in `vfa_example_script.m`

References: 

```
Zhao, L., Mulkern, R., Tseng, C.-H., Williamson, D., Patz, S., Kraft, R., Walsworth, R. L., Jolesz, F. A., and Albert, M. S. Gradient-echo imaging considerations for hyperpolarized 129xe mr. J Magn Reson B 113, 2 (November 1996), 179–183.
```
```
Xing Y, Reed GD, Pauly JM, Kerr AB, and Larson PEZ. Optimal Variable Flip Angle Schemes For Dynamic Acquisition Of Exchanging Hyperpolarized Substrates.  J Magn Reson. 2013 Sep; 234:75-81.  http://www.ncbi.nlm.nih.gov/pubmed/23845910
```
### optimal_SNR_angles()

Optimization of total metabolic product signal when summed across all acquisitions.
Results are nearly identical to `vfa_opt_signal`, as demonstrated in `optimal_SNR_demo.m`
