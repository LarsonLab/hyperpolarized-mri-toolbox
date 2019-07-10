function kspace_all = recon_cs3d_dyn_sivic_v1(root_dir, fb_root_name, samp_pattern_name, ...
    numReps, Navg, weights, ItParams)

% function to do an L1 reconstruction of 3D dynamic data
% currently using wavelet in time as sparsifying transform
%
% PEZL 2.9.09, 2.9.12

addpath /netopt/share/lib/local/brain/matlab
addpath(genpath('/home/plarson/matlab/reconstruction/rand_epsi/'))
addpath(genpath('/home/plarson/matlab/reconstruction/sparseMRI/'))
addpath(genpath('/home/plarson/matlab/utilities/Wavelab850'))

if nargin < 5
  Navg = 1;
end

if nargin < 6  
  weights = [1e-6 1e-2];
end

if nargin < 7
  ItParams = [16 6];
end

plot_test = 0;


load([root_dir samp_pattern_name]); % File from blip generation

S = size(loc_samp_3d_dyn);
length_x = S(1); length_y = S(2); length_f = S(3);
length_z = 16;

% data_all = zeros(length_f, length_x, length_y, length_z, numReps);
mask_all = zeros(length_f, length_x, length_y, numReps);


for n = 1:numReps
    Isum = [1:Navg] + (n-1)*Navg;
    
    mask_all(:,:,:,n) = ...
        shiftdim(any(loc_samp_3d_dyn(:,:,:,Isum),4),2);
end

load([root_dir fb_root_name]); % File from 3d undersample dataset

timedim = 4; % after squeeze

%%%%%%%%%%%%%
% L1 Recon
%%%%%%%%%%%%%
Itnlim = ItParams(1);
loopIter = ItParams(2);
TVWeight = weights(1); %1e-3; %1e-2;
xfmWeight = weights(2); %1e-2; %2e-2;;
wavelet_time_scale = 2;  % similar between 0,1,2

N = size(mask_all); % k-space and image size do not have to be the same
NN1 = 2^ceil(log(N(1))/log(2));
if numReps > 1
    NN4 = 2^ceil(log(N(4))/log(2)); % interpolate to diadic number for wavelets
else
    NN4 = 1;
end
NN = [NN1,N(2),N(3),NN4];

data_hybrid = fftshift(ifft(data_all, [], 4), 4);
kspace_hybrid = zeros(length_f, length_x, length_y, length_z, numReps);

poolobj = parpool;

parfor z = 1:length_z
    disp(sprintf('====== START: slice %d =======', z))

    XFM = Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);  % wavelet along time - seems to be the most sparse
    FT = pnDFT_time(mask_all, NN, timedim, 1, 2);

    param = init;
    param.FT = FT;
    param.XFM = XFM;
    param.TV = TVOP_3D;  %or TVOP_4D?  not much difference
    param.Itnlim = Itnlim;
    param.TVWeight = TVWeight;  % TV penalty
    param.xfmWeight = xfmWeight;  % L1 wavelet penalty

    data_slice = squeeze(data_hybrid(:,:,:,z,:));
    im_init = FT'*data_slice;

    % scale k-space data
    scale = max(abs(im_init(:)));
    param.data = (1/scale) * data_slice;  % Fourier data
    
    % input should be in sparse domain
    im_init = FT'*(param.data);
    resrecon = XFM*im_init;
%    res_init = resrecon;

    tic
    for n = 1:loopIter
        % enfore data consistency
        data_res = crop(fftnc_time(XFM' * resrecon, timedim), N);
        data_res = param.data.*mask_all + data_res.*(~mask_all);
        resrecon = XFM * ifftnc_time(zpad(data_res, NN), timedim);

        resrecon = fnlCg(resrecon,param);
    end

    toc

    srecon = (XFM' * resrecon);
    Srecon = fftnc_time(srecon, timedim);
    kspace_hybrid(:,:,:,z,:) = crop(Srecon*scale, N);

    disp(sprintf('z = %d\n',z));

end

kspace_all = fft(fftshift(kspace_hybrid, 4), [], 4);

delete(poolobj)

end


