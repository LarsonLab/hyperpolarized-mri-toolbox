% fit slab dynamic MRS data from repeated injections ("shots") in a single rat

% Data and studies described in
% Hu, Simon, Peder E Z Larson, Mark Vancriekinge, Andrew M Leach, Ilwoo Park, Christine Leon,  et al.
% "Rapid Sequential Injections of Hyperpolarized [1-13C]pyruvate in Vivo
% Using a Sub-Kelvin, Multi-Sample DNP Polarizer"
% Magn Reson Imaging 31, no. 4 (May 2013): 490-96.
% https://doi.org/10.1016/j.mri.2012.09.002.

clear all
plot_flag = 1;
int_range = 25;
R1P_est = 1/30; R1L_est = 1/25;  R1B_est = 1/15; R1A_est = 1/25;
kPL_est = 0.01; kPB_est = .01; kPA_est = .01;

for ratnum = 1:5
    
    switch ratnum
        case 1
            Nshots = 4;
            lac_center = 970;
            pyr_hyd_center = 1020;
            ala_center = 1056;
            pyr_center = 1130;
            bicarb_center = 1262;
        case 2
            Nshots = 3;
            lac_center = 967;
            pyr_hyd_center = 1018;
            ala_center = 1053;
            pyr_center = 1128;
            bicarb_center = 1260;
        case 3
            Nshots = 3;
            lac_center = 968;
            pyr_hyd_center = 1019;
            ala_center = 1054;
            pyr_center = 1129;
            bicarb_center = 1260;
        case 4
            Nshots = 3;
            
            lac_center = 876;
            pyr_hyd_center = 927;
            ala_center = 963;
            pyr_center = 1037;
            bicarb_center = 1168;
            
        case 5
            Nshots = 3;
            lac_center = 875;
            pyr_hyd_center = 926;
            ala_center = 961;
            pyr_center = 1036;
            bicarb_center = 1168;
    end
    
    
    clear params_fixed params_est params_fit kPL kPB kPA
    params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
    params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
    params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;
    
    for n = 1:Nshots
        load(sprintf('rat%d/shot%d',ratnum,n));
        pyr = real(sum(spectra_dynamic(pyr_center+[-int_range:int_range],:),1));
        lac = real(sum(spectra_dynamic(lac_center+[-int_range:int_range],:),1));
        ala = real(sum(spectra_dynamic(ala_center+[-int_range:int_range],:),1));
        bicarb = real(sum(spectra_dynamic(bicarb_center+[-int_range:int_range],:),1));
        
        % only fitting kPL - gives similar results
        %         [params_fit Sfit] = fit_pyr_kinetics([pyr;lac] , TR, repmat(flip, [2 length(pyr)]), params_fixed, params_est, [], plot_flag);
        %         kPL_only(n) = params_fit.kPL;
        
        [params_fit Sfit] = fit_pyr_kinetics([pyr;lac;bicarb;ala] , TR, repmat(flip, [4 length(pyr)]), params_fixed, params_est, [], plot_flag);
        kPL(n) = params_fit.kPL; kPB(n) = params_fit.kPB; kPA(n) = params_fit.kPA;
        
    end
    
    disp(['Rat #' int2str(ratnum)]);
    disp(['   kPL fits: ' num2str(kPL,3) ' 1/s'])
    disp(['   kPB fits: ' num2str(kPB,3) ' 1/s'])
    disp(['   kPA fits: ' num2str(kPA,3) ' 1/s'])
end