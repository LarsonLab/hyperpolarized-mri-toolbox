function [kTRANS, kMaps_out, Mz0Maps_out, metImages, w] = brainweb_metabolic_phantom(kineticRates, ktransScales, Mz0, sampSize, outputSize, simParams, inputFunction, isFuzzy, linear_kTRANS_grad, augmentParams, brain_idx, augmentSeed)
% BRAINWEB_METABOLIC_PHANTOM generates standardized 3-dimensional perfusion
%   and metabolism maps for simulated experiments. Supports 3 chemical pool
%   kinetic rate mapping.
%
%   Toolboxes required: Image Processing
%
%   Parameters:
%       kineticRates    = the kinetic rates to simulate, [(# of chemical
%                         pools) - 1, 3 (tissue types: vasc, gm, wm)]
%       ktransScales    = volume transfer constants for different tissue
%                         compartments [vasculature, gray matter, white matter]
%                         if linear_kTRANS_grad=true, should include low
%                         and high kTRANS [[vasc_low, GM_low, WM_low]; [vasc_hi, GM_hi, WM_hi]]
%       Mz0             = Initial magnetization per metabolite per compartment,
%                         size [Nmets 3] in order vessels, gm, wm
%       sampSize        = 3x3 or 1x3 vector for desired sampling ("acquisition") matrix size of each
%                         dimension for all mets
%                         [met dim(nx,ny,nz)], default = [16 16 8]
%       outputSize      = 1x3 vector for desored output matrix size,
%                         default = [64 64 8]
%       simParams       = parameters used for the kinetic simulations: 
%                         Tarrival, Tbolus, TR, Nt, R1, flips; if empty won't
%                         generate metImages, default = empty struct
%       inputFunction  = [1 Nt] vector as input into each voxel
%
%  Optional Parameters:
%       isFuzzy         = logical flag to indicate if fuzzy tissue boundaries
%                         will be used, default = true
%  linear_kTRANS_grad   = boolean flag for whether kTRANS should be set to 
%                         a gradient
%       augmentParams   = random augmentation parameters XTranslation, 
%                         YTranslation etc, if not defined a single 
%                         unaugmented phantom will be generated,
%                         default = empty struct
%       brain_idx       = specify if want to use a specific brain (1-20),
%                         default = randomly choose a brain
%       augmentSeed     = random seed for augmentations, default no seed
%
%   Outputs:
%       kTRANS      = generated perfusion map
%       kMaps       = generated rate maps for 1->2 and 1->3
%       Mz0Maps     = generated Mz0 maps
%       metImages   = simulated metabolite dynamic images
%       w           = simulated coil sensitivity maps
%
%   Author:
%       Jasmine Hu
%       Anna Bennett
%       Sule Sahin
%
% Copyright, 2024

    % parse input arguments
    arguments
        kineticRates (:,3) double {mustBeNumeric} = [0.1, 0.2, 0.3; 0, 0, 0]
        ktransScales (:,3) double {mustBeNumeric} = [1, 0.3, 0.3]
        Mz0 (:,3) double {mustBeNumeric} = [0, 0, 0]
        sampSize (:,3) double {mustBeInteger} = [16, 16, 8]
        outputSize (1,3) double {mustBeInteger} = [64, 64, 8]
        simParams struct = struct([])
        inputFunction double = []
        isFuzzy double {mustBeNumericOrLogical} = true
        linear_kTRANS_grad double {mustBeNumericOrLogical} = false
        augmentParams struct = struct()
        brain_idx {mustBeInteger} = [] %TODO: should be int between 1-19
        augmentSeed double {mustBeInteger, mustBePositive, mustBeNonzero} = []
    end

    if augmentSeed
        rng(augmentSeed)
    else
        rng("shuffle")
    end

    if and(linear_kTRANS_grad, size(ktransScales,1)==1)
        ktransScales = repmat(ktransScales, [2 1]);
        warning("kTRANS linear gradient is 'true' but low/high kTRANS values not defined. Using user-defined kTRANS values as both low and high.")
    elseif and(~linear_kTRANS_grad, size(ktransScales,1)>1)
        ktransScales = squeeze(ktransScales(1,:));
        warning("kTRANS linear gradient is 'false' but low/high kTRANS values defined. Using first row of kTRANS values.")
    end

    nMets = size(kineticRates,1) + 1;

    if size(sampSize) == [1 3]
        maxSampSize = sampSize;
        sampSize = repmat(sampSize, [nMets 1]);
    else
        [~,i] = max(sampSize(:,1));
        maxSampSize = sampSize(i,:);
    end

    % add resources dir to path
    fileDir = split(mfilename('fullpath'),'/');
    utilDir = fullfile(string(join(fileDir(1:end-1),'/')),'/util');
    addpath(utilDir)
    
    % define default augmentation parameters and set augmentParams to them
    % if not defined
    defaultAugParams.XReflection = false; defaultAugParams.YReflection = false;
    defaultAugParams.Rotation = [0 0]; defaultAugParams.Scale = [1 1];
    defaultAugParams.XShear = [0 0]; defaultAugParams.YShear = [0 0];
    defaultAugParams.XTranslation = [0 0]; defaultAugParams.YTranslation = [0 0];
    augs = fields(defaultAugParams);
    for f=1:size(augs,1)
        if ~isfield(augmentParams,augs{f})
            augmentParams(1).(augs{f}) = defaultAugParams.(augs{f});
        end
    end

    % load base anatomical information
    if isempty(brain_idx)
        brain_idx = randi(19);
    end
    if isFuzzy
        brainwebFile = 'brainweb_fuzzy.mat';
    else
        brainwebFile = 'brainweb.mat';
    end
    baseMaskFile = dir(fullfile(utilDir,num2str(brain_idx),brainwebFile));
    if isempty(baseMaskFile)
        error('Error. \nBrainWeb mask file, %s, not found in resources directory.',brainwebFile)
    end
    
    load(fullfile(baseMaskFile.folder,brainwebFile),'im_mask');
    vasc_mask = squeeze(im_mask(:,:,:,1));
    gm_mask = squeeze(im_mask(:,:,:,2));
    wm_mask = squeeze(im_mask(:,:,:,3));
    brain_mask = vasc_mask + gm_mask + wm_mask;
    maskSize = size(im_mask,1:3);
    
    weights = coil_dist_map(brain_mask, simParams.coil_lim);

    nTissues = 3;
    if size(im_mask,4) ~= nTissues
        error('Unexpected number of tissues present in the imported masks, ask for help, idk.');
    end
    
    % parameters for output map generation
    permuted_mask = double(permute(im_mask,[4 1 2 3]));
    sumWeights = sum(im_mask,4);

    if linear_kTRANS_grad 
        % generate kTRANS gradients per compartment
        grad_vasc = generate_linear_gradient(maskSize, ktransScales(1,1), ktransScales(2,1));
        grad_gm = generate_linear_gradient(maskSize, ktransScales(1,2), ktransScales(2,2));
        grad_wm = generate_linear_gradient(maskSize, ktransScales(1,3), ktransScales(2,3));
        
        % generate the kTRANS masked volume
        kTRANS_vasc = squeeze(permuted_mask(1,:,:,:)) .* grad_vasc;
        kTRANS_gm = squeeze(permuted_mask(2,:,:,:)) .* grad_gm;
        kTRANS_wm = squeeze(permuted_mask(3,:,:,:)) .* grad_wm;
        kTRANS_wSum = kTRANS_vasc + kTRANS_gm + kTRANS_wm;
        kTRANS = squeeze(kTRANS_wSum)./sumWeights;
        kTRANS(isnan(kTRANS)) = 0;
    else
        kTRANS = create_map(permuted_mask, ktransScales, sumWeights);
    end
    
    % generate the kinetic rate maps
    k_1_2_MAP = create_map(permuted_mask, kineticRates(1,:), sumWeights);
    k_1_3_MAP = create_map(permuted_mask, kineticRates(2,:), sumWeights);

    % generate Mz0 maps
    Mz0P_MAP = create_map(permuted_mask, Mz0(1,:), sumWeights);
    Mz0L_MAP = create_map(permuted_mask, Mz0(2,:), sumWeights);
    Mz0B_MAP = create_map(permuted_mask, Mz0(3,:), sumWeights);

    %adjust FOV in z
    cropidx1 = randi([40 60]); % TODO: make this optional?
    cropidx2 = randi([300 320]);
    kTRANS = kTRANS(:,:,cropidx1:cropidx2);
    k_1_2_MAP = k_1_2_MAP(:,:,cropidx1:cropidx2);
    k_1_3_MAP = k_1_3_MAP(:,:,cropidx1:cropidx2);
    Mz0P_MAP = Mz0P_MAP(:,:,cropidx1:cropidx2);
    Mz0L_MAP = Mz0L_MAP(:,:,cropidx1:cropidx2);
    Mz0B_MAP = Mz0B_MAP(:,:,cropidx1:cropidx2);
    w = weights(:,:,cropidx1:cropidx2);
    
    % resample/downsample maps to desired SAMPLE size
    kTRANS = imresize3(kTRANS, maxSampSize);
    k_1_2_MAP = imresize3(k_1_2_MAP, maxSampSize);
    k_1_3_MAP = imresize3(k_1_3_MAP, maxSampSize);
    Mz0P_MAP = imresize3(Mz0P_MAP, maxSampSize);
    Mz0L_MAP = imresize3(Mz0L_MAP, maxSampSize);
    Mz0B_MAP = imresize3(Mz0B_MAP, maxSampSize);
    w = imresize3(w, maxSampSize);
   
    % Random Augmentations (requires image processing toolbox)
    tform = randomAffine2d('Rotation',augmentParams.Rotation, 'Scale', ...
        augmentParams.Scale, 'XReflection', augmentParams.XReflection, ...
        'YReflection', augmentParams.YReflection, 'XTranslation', ...
        augmentParams.XTranslation, 'YTranslation', augmentParams.YTranslation, ...
        'XShear', augmentParams.XShear, 'YShear', augmentParams.YShear);
    outputView = affineOutputView(size(kTRANS),tform);
    kTRANS = imwarp(kTRANS,tform,OutputView=outputView); 
    k_1_2_MAP = imwarp(k_1_2_MAP,tform,OutputView=outputView);
    k_1_3_MAP = imwarp(k_1_3_MAP,tform,OutputView=outputView);
    Mz0P_MAP = imwarp(Mz0P_MAP,tform,OutputView=outputView);
    Mz0L_MAP = imwarp(Mz0L_MAP,tform,OutputView=outputView);
    Mz0B_MAP = imwarp(Mz0B_MAP,tform,OutputView=outputView);
    w = imwarp(w,tform,OutputView=outputView);
    
    %flip 3rd dimension to match in vivo convention
    kTRANS = flip(kTRANS,3);
    kMaps = flip(cat(4,k_1_2_MAP,k_1_3_MAP),3); 
    Mz0Maps = flip(cat(4,Mz0P_MAP,Mz0L_MAP, Mz0B_MAP),3);
    w = flip(w,3);
    w = w ./ max(w(:)); % normalize coil sens weights to max of 1 after imresize

    if ~isempty(simParams) % output simulated metabolite dynamic images
        % simulate signals
        % store simulation parameters a in a struct,
        % eventually this should be a custom class

        if isempty(inputFunction)
            inputFunction = zeros([1 simParams.Nt]);
        end
        
        % first generate voxelwise dynamics
        metImages_sp = zeros(cat(2,maxSampSize,[nMets, simParams.Nt]));
        metImages = zeros(cat(2,outputSize,[nMets, simParams.Nt]));
        for Ix = 1:maxSampSize(1)
            for Iy = 1:maxSampSize(2)
                for Iz = 1:maxSampSize(3)
                    Mz0_vx = [Mz0Maps(Ix, Iy, Iz, 1) Mz0Maps(Ix, Iy, Iz, 2) Mz0Maps(Ix, Iy, Iz, 3)];
                    [metImages_sp(Ix,Iy,Iz,:,:), ~] = simulate_Nsite_model(Mz0_vx, simParams.R1, [kMaps(Ix,Iy,Iz,1) 0; kMaps(Ix,Iy,Iz,2) 0], simParams.flips, simParams.TR, inputFunction*kTRANS(Ix,Iy,Iz) );
                end
            end
        end
        
        % multiply by coil sens weights
        metImages_sp = metImages_sp .* repmat(w, [1 1 1 size(metImages_sp,4) size(metImages_sp,5)]); 
        
        % multi res capability
        metImages_mres = cell(nMets,1); % cast metImages into cell to allow different matrix sizes
        for Imet=1:nMets

            % add rician noise
            std_noise = max(sum(squeeze(metImages_sp(:,:,:,Imet,:)),4),[],'all') ./ (simParams.SNR(Imet) * sqrt(simParams.Nt));
            noise_R = randn(cat(2,sampSize(Imet,:),[simParams.Nt]))* std_noise; 
            noise_I = randn(cat(2,sampSize(Imet,:),[simParams.Nt]))* std_noise;
            metImages_w_noise = sqrt((metImages(nmet)+ noise_R).^2 + noise_I.^2);

            if maxSampSize == sampSize(Imet,:)
                temp = squeeze(metImages_sp(:,:,:,Imet,:));
                metImages_mres{Imet} = sqrt((temp + noise_R).^2 + noise_I.^2);
            else
                temp=zeros(cat(2,sampSize(Imet,:),simParams.Nt));
                for t=1:simParams.Nt
                    temp(:,:,:,t) = imresize3(metImages_sp(:,:,:,Imet,t), sampSize(Imet,:), 'box'); %box downsampling will effectively average across voxels
                end
                metImages_mres{Imet} = sqrt((temp + noise_R).^2 + noise_I.^2);
                clear temp;
            end
        end

        % resize to desired OUTPUT size
        for Imet = 1:size(metImages_sp,4)
            temp = metImages_mres{Imet};
            for It = 1:size(metImages_sp,5)
                metImages(:,:,:,Imet,It) = imresize3(squeeze(temp(:,:,:,It)), outputSize, 'lanczos3');
                %metImages(:,:,:,Imet,It) = zeropad(squeeze(temp(:,:,:,It)), outputSize(1:2));
            end
        end
    
    else
        metImages = 0;
    end

    % resample maps to desired OUTPUT size
    kTRANS = imresize3(kTRANS, outputSize, 'lanczos3');
    kMaps_out = zeros(cat(2,outputSize,nMets-1));
    Mz0Maps_out = zeros(cat(2,outputSize,nMets));
    for n=1:nMets-1
        kMaps_out(:,:,:,n) = imresize3(kMaps(:,:,:,n), outputSize, 'lanczos3');
    end
    for n=1:nMets
        Mz0Maps_out(:,:,:,n) = imresize3(Mz0Maps(:,:,:,n), outputSize, 'lanczos3');
    end

end

function [grad] = generate_linear_gradient(maskSize, kTRANS_low, kTRANS_high)

    x = linspace(-1, 1, maskSize(1));
    y = linspace(-1, 1, maskSize(2));
    z = linspace(-1, 1, maskSize(3));
    [X, Y, Z] = meshgrid(x, y, z);
    grad = 0.5*(kTRANS_high - kTRANS_low)*Y + 0.5*(kTRANS_low + kTRANS_high);
end

function [map] = create_map(mask, rates, sumWeights)
    map_wSum = pagemtimes(rates,mask);
    map = squeeze(map_wSum)./sumWeights;
    map(isnan(map)) = 0;
end