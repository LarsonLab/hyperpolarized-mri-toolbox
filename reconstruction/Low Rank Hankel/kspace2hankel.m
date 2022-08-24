function [onehankelmtx_mcoil,hankel_params] = kspace2hankel(kspace_hybrid,forwardmtx,hankel_params)
% construct Hankel matrix from 3D EPSI hybrid k-space data (k-space + dyn)
% v1. 20171018

onehankelmtx_mcoil = [];

numReps = hankel_params.numReps;
numCoils = hankel_params.numCoils;

for cc = 1:numCoils
    data = kspace_hybrid(:,:,:,:,cc);
    dim_data = size(data);
    hankelmtx = data(forwardmtx);

    % concatenate along dynamic dimension
    hankelmtx_temp1 = [];
    for i_dyn = 1:numReps,
        hankelmtx_temp1 = cat(2,hankelmtx_temp1, squeeze(hankelmtx(:,:,i_dyn)));
    end
    onehankelmtx_mcoil = cat(1,onehankelmtx_mcoil,hankelmtx_temp1);  %cat along coil dim
end

hankel_params.dim_data = dim_data;
hankel_params.dim_hankel = size(forwardmtx);

% varargout{1} = hankel_params;

end