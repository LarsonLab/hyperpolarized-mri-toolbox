function kspace_all = recon_cs3d_for_dyn_linearized_Bregman_dynsvd_slidekxy_mcoil(root_dir, fb_root_name, samp_pattern_name, ...
    numReps, numCoils, Navg, weights, ItParams)

% function to do an L1 reconstruction of 3D dynamic data
% currently using wavelet in time as sparsifying transform
%
% PEZL 2.9.09, 2.9.12
% pcao 07.03.16, multi coil recon

addpath /netopt/share/lib/local/brain/matlab
addpath /home/pcao/matlab/utilities/rsvd
addpath /home/pcao/matlab/utilities/lmsvd
addpath(genpath('/home/plarson/matlab/reconstruction/rand_epsi/'))
addpath(genpath('/home/plarson/matlab/reconstruction/sparseMRI/'))
addpath(genpath('/home/plarson/matlab/utilities/Wavelab850'))
addpath(genpath('/home/pcao/matlab/utilities/svdandpca'))
cd(root_dir)
if nargin < 6
  Navg = 1;
end

if nargin < 7  
  weights = 300;  %threshold on svd
end

if nargin < 8
  ItParams = [300 150]; %num of singular values, num of iterations
end

wldyn = 1;


load([root_dir samp_pattern_name]); % File from blip generation

S = size(loc_samp_3d_dyn);

length_x = S(1); length_y = S(2); length_f = S(3);

% data_all = zeros(length_f, length_x, length_y, length_z, numReps, numCoils);
% mask_all = zeros(length_f, length_x, length_y, numReps, numCoils);

for cc = 1:numCoils
    if numCoils > 1
        coilnum = sprintf('_%d',cc); %if multichannel coil the names are usually ssssss_1.ddf, ssssss_2.ddf.....
    else
        coilnum = [];
    end
    for n = 1:numReps
        dataRep{n} = read_ddf_image([sprintf('%s%s%02d', root_dir , fb_root_name , n) coilnum]);
        
        Isum = [1:Navg] + (n-1)*Navg;
        mask_all(:,:,:,n,cc) = ...
            shiftdim(any(loc_samp_3d_dyn(:,:,:,Isum),4),2);
        data_all(:,:,:,:,n,cc) = ...
            dataRep{n}.img(:,:,:,:);    
    end
end


length_z = size(data_all,4);
timedim = 4; % after squeeze

%%%%%%%%%%%%%
% svd Recon
%%%%%%%%%%%%%
loopIter = ItParams(2);
wavelet_time_scale = 2;  % similar between 0,1,2

N = size(mask_all); % k-space and image size do not have to be the same

NN1 = 2^ceil(log(N(1))/log(2));
if numReps > 1
    NN4 = 2^ceil(log(N(4))/log(2)); % interpolate to diadic number for wavelets
else
    NN4 = 1;
end

if wldyn == 0
NN4 = numReps;
end

NN = [NN1,N(2),N(3),NN4,N(5)];

data_hybrid = fftshift(ifft(data_all, [], 4), 4);%fft along slice dim
kspace_hybrid = zeros(length_f, length_x, length_y, length_z, numReps, numCoils);


% poolobj = parpool;

%generate the LUT for forward and backward matrix
%truncate of data to accelerate the svd
t_trunc = NN4;
sp_trunc = NN1;

sp_trunc_range = (NN1 - sp_trunc)/2 + 1:(NN1 + sp_trunc)/2 ;
list = exist([root_dir '/LUT_forward_backward.mat'], 'file');
if ~list
    forwardmtx = gen_forwardmtx_slidexy_col(zeros([sp_trunc,N(2),N(3),t_trunc]),10);% sp_trunc-16
    [backwardmtx, sumdensity] = gen_backwardmtx2(zeros([sp_trunc,N(2),N(3),t_trunc]), forwardmtx);
    save([root_dir '/LUT_forward_backward.mat'], 'backwardmtx', 'sumdensity', 'forwardmtx', '-v7.3');
else
    load([root_dir '/LUT_forward_backward.mat']);
end

% scale k-space data
scale = max(abs(data_hybrid(:)));
    
% parfor z = 1:length_z
for z = 1:length_z
    disp(sprintf('====== START: slice %d =======', z))
    
    if wldyn == 1
        XFM = Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);  % wavelet along time - seems to be the most sparse
    else
        XFM = 1;
    end
     
    % scale k-space data
    org_data = (1/scale) * reshape(data_hybrid(:,:,:,z,:,:),N);  % Fourier data

    data_res = zeros(size(org_data));
    data_res_ex = zeros(size(org_data));
    srecon = zeros(NN);
    tic
    for n = 1:loopIter
        pre_data_res_ex = data_res_ex; %for kicking
        data_res_ex = crop(srecon, N);    
        data_res = data_res + (org_data - data_res_ex).*abs(mask_all);
        resrecon = XFM * zpad(data_res, NN);   %cp:wavelet, sparse transform

        spec = sum(sum(sum(sum(resrecon,2),3),4),5);
        spec_recon = sum(sum(sum(sum(pre_data_res_ex,2),3),4),5);

%         figure(911)
%         subplot(1,2,1), plot(real(spec));
%         subplot(1,2,2), plot(real(spec_recon));
%         drawnow;
         
        [resrecon dim_data dim_hankel rank] = SVD_threshold_hankel_dynsvd_slidekxy_dynave_mcoil_v2(resrecon,'fixed_th',ItParams(1), forwardmtx, backwardmtx, sumdensity,weights, numReps, numCoils); %working version
        
        res = abs(org_data-data_res_ex).*mask_all;
        norm_res = sum(res(:).^2)^1/2;
        absdata = abs(org_data).*mask_all;
        norm_data = sum(absdata(:).^2)^1/2;
        disp(sprintf('slice = %d; loopIter = %d; MSE_residual = %d', ...
            z, n, norm_res/norm_data));
        
        srecon = (XFM' * resrecon);  %cp:srecon  kx    ky    t
            
        if (norm_res/norm_data < 10e-2)
            break;
        end

    end

    toc
    
    kspace_hybrid(:,:,:,z,:,:) = crop(srecon*scale, N);
    imag_hybrid(:,:,:,z,:,:) = crop(XFM' * resrecon*scale, N);

%     disp(sprintf('z = %d\n',z));

end

% save image data
kspace_all = fft(fftshift(kspace_hybrid, 4), [], 4);

% create new ddf and cmplx
for cc = 1:numCoils
    if numCoils > 1
        coilnum = sprintf('_%d',cc); %if multichannel coil the names are usually ssssss_1.ddf, ssssss_2.ddf.....
    else
        coilnum = [];
    end
    for n = 1:numReps
        dataRep{n}.img = kspace_all(:,:,:,:,n,cc);
        write_ddf_image([sprintf('%s%s%02d', root_dir, fb_root_name, n) coilnum],dataRep{n});
    end
end

end
