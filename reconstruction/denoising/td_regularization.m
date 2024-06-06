function imgOut = td_regularization(imgIn,tensorRank,plot_flag)
% based on the following paper
% "PET by MRI: Glucose Imaging by 13C-MRS without Dynamic Nuclear 
% Polarization by Noise Suppression through Tensor Decomposition 
% Rank Reduction", Brender, JR et al., bioRxiv 2018

% imgIn: input n-dimensional image array, e.g. imgIn(f,x,y,z,dyn,...etc)
% imgOut: denoised n-dimensional image array
% tensorRank: n-element vector that specifies # of ranks to keep in each
%             spectral or spatial dimension.
%             (similar to eigenvalues in conventional SVD)
%
% Note this function requires the tensorlab toolbox in
% path (https://www.tensorlab.net/index.html)
%
% v1: 20180822 HYC

if nargin < 3 || isempty(plot_flag)
    plot_flag = 1;
end

[U,S,sv] = mlsvd(imgIn);

ndim_imgIn = ndims(imgIn);
size_S_lowrank_temp = size(imgIn);
S_lowrank = S;
% tensor singular value thresholding
for i_dim = 1:ndim_imgIn
    U_lowrank{i_dim} = U{i_dim}(:, 1:tensorRank(i_dim));
%     S_lowrank_1 = S(1:tensorRank(1), 1:tensorRank(2), 1:tensorRank(3), 1:tensorRank(4));
    S_lowrank = S_lowrank(1:tensorRank(i_dim),:);
    S_lowrank = reshape(S_lowrank,[tensorRank(i_dim) size_S_lowrank_temp(2:end)]);
    size_S_lowrank_temp = [size_S_lowrank_temp(2:end) tensorRank(i_dim)];
    S_lowrank = shiftdim(S_lowrank,1);
end
imgOut = lmlragen(U_lowrank, S_lowrank);

% plot tensor values in each dimension
if plot_flag == 1
    figure,
    for i_dim = 1:ndim_imgIn
      y = sv{i_dim};
      subplot(1,ndim_imgIn,i_dim);
      semilogy(y);
      title(sprintf('Dim %01d',i_dim));
    end
end

end