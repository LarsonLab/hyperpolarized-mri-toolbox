function imgOut = td_regularization_fast(imgIn_temp,tensorRank,plot_flag)
% based on the following paper
% "PET by MRI: Glucose Imaging by 13C-MRS without Dynamic Nuclear 
% Polarization by Noise Suppression through Tensor Decomposition 
% Rank Reduction", Brender, JR et al., bioRxiv 2018

% imgIn: input n-dimensional image array, e.g. imgIn(f,x,y,z,dyn,...etc)
% imgOut: denoised n-dimensional image array
% tensorRank: n-element vector that specifies # of ranks to keep in each
%             spectral or spatial dimension.
%             (similar to eigenvalues in conventional SVD)

% v1: 20180822 HYC

if nargin < 3 || isempty(plot_flag)
    plot_flag = 1;
end

U = imgIn_temp.U;
S = imgIn_temp.S;
sv = imgIn_temp.sv;
ndim_imgIn = imgIn_temp.ndim_imgIn;

range_tensorRank = cell(1,ndim_imgIn);
% tensor singular value thresholding
for i_dim = 1:ndim_imgIn
    U_lowrank{i_dim} = U{i_dim}(:, 1:tensorRank(i_dim));
    range_tensorRank{i_dim} = 1:tensorRank(i_dim);
end
% S_lowrank = S(1:tensorRank(1), 1:tensorRank(2), ... , 1:tensorRank(n));
S_lowrank = S(range_tensorRank{:});
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