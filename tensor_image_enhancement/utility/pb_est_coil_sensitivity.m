function [svdWeights,varargout] = pb_est_coil_sensitivity(Data_5D_orig)
    % params
    size_Data_5D_orig = size(Data_5D_orig); % [f,x,y,coil,dyn]
    Data_5D_orig_AUC = sum(Data_5D_orig,5);
    noiseMask = logical([zeros(225,1); ones(31,1)]);
    
    svd_spectra = zeros(size_Data_5D_orig(1),size_Data_5D_orig(2),size_Data_5D_orig(3));
    svdQuality = zeros(size_Data_5D_orig(2),size_Data_5D_orig(3));
    svdCoilAmpt = zeros(size_Data_5D_orig(4),size_Data_5D_orig(2),size_Data_5D_orig(3));
    svdWeights = zeros(size_Data_5D_orig(4),size_Data_5D_orig(2),size_Data_5D_orig(3));
    % cal sensitivities using SVD
    for x = 1:size_Data_5D_orig(2),
    for y = 1:size_Data_5D_orig(3),
        [svdCoilAmpt(:,x,y),svdWeights(:,x,y),svd_spectra(:,x,y), svdQuality(x,y)] = ...
            svdRecon(squeeze(Data_5D_orig_AUC(:,x,y,:)), noiseMask, 'noiseCov', 'disable');  
    end
    end
    varargout{1} = svdCoilAmpt;
end