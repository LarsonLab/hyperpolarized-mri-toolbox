function [kspace_hybrid] = hankel2kspace(onehankelmtx_mcoil,backwardmtx,hankel_params)
% extract 3D EPSI hybrid k-space data from Hankel matrix
% v1. 20171018

dim_data = hankel_params.dim_data;
dim_hankel = hankel_params.dim_hankel;
numReps = hankel_params.numReps;
numCoils = hankel_params.numCoils;
sumdensity = hankel_params.sumdensity;

kspace_hybrid = zeros([dim_data numCoils]);

for cc = 1:numCoils,
    hankelmtx = zeros(dim_hankel);%
    for i_dyn = 1:numReps,
        coil_dim_idx = (1:dim_hankel(1)) + (cc-1) * dim_hankel(1);%extract along coil dim
        dyn_dim_idx = (1:dim_hankel(2)) + (i_dyn - 1) * dim_hankel(2);
        hankelmtx(:,:,i_dyn) = onehankelmtx_mcoil(coil_dim_idx,dyn_dim_idx);
    end
    hankelmtx(1) = 0;
    backwardmtx(backwardmtx == 0) = 1;
    realinemtx = hankelmtx(backwardmtx);
    kspace_hybrid(:,:,:,:,cc) = reshape( sum(realinemtx, 2)./sumdensity, dim_data);
end

end