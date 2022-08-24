function [res dim_data dim_hankel rank diagv] = SVD_threshold_hankel_dynsvd_slidekxy_dynave_mcoil_v2(datain,mode,k, forwardmtx, backwardmtx, sumdensity,th, numReps, numCoils)
%fast hankel matrix projection for SVD threshold method

ss = nnz(abs(datain(:)));

if ss < k
    k = 1;
end
hankel_params.numReps = numReps;
hankel_params.numCoils = numCoils;
hankel_params.sumdensity = sumdensity;
% k space to hankel
[onehankelmtx_mcoil,hankel_params] = kspace2hankel(datain,forwardmtx,hankel_params);
% SVD threshold
[onehankelmtx_mcoil(:,:), rank,diagv] = SVD_thres_kxy( onehankelmtx_mcoil(:,:), mode, k, th);
% hankel to k space
res = hankel2kspace(onehankelmtx_mcoil,backwardmtx,hankel_params);

dim_data = hankel_params.dim_data;
dim_hankel = hankel_params.dim_hankel;
 
end

