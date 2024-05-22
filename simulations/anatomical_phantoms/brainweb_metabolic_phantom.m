function [kTRANS, kMaps, Mz0Maps, metImages] = brainweb_metabolic_phantom(kineticRates, ktransScales, Mz0, matSize, simParams, inputFunction, isFuzzy, linear_kTRANS_grad, augmentParams, brain_idx, augmentSeed)
% BRAINWEB_METABOLIC_PHANTOM generates standardized 3-dimensional perfusion
%   and metabolism maps for simulated experiments. Supports 3 chemical pool
%   kinetic rate mapping.
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
%       matSize         = 1x3 vector for desired matrix size of each dimension
%                         [nx ny nz], default = [16 16 8]
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
        matSize (1,3) double {mustBeInteger} = [16, 16, 8]
        simParams struct = struct([])
        inputFunction double = []
        isFuzzy double {mustBeNumericOrLogical} = true
        linear_kTRANS_grad double {mustBeNumericOrLogical} = false
        augmentParams struct = struct()
        brain_idx {mustBeInteger} = []
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
        brain_idx = randi(20);
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
    
    weights = coil_dist_map(brain_mask);

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
    cropidx1 = randi([40 60]); % make this optional?
    cropidx2 = randi([300 320]);
    kTRANS = kTRANS(:,:,cropidx1:cropidx2);
    k_1_2_MAP = k_1_2_MAP(:,:,cropidx1:cropidx2);
    k_1_3_MAP = k_1_3_MAP(:,:,cropidx1:cropidx2);
    Mz0P_MAP = Mz0P_MAP(:,:,cropidx1:cropidx2);
    Mz0L_MAP = Mz0L_MAP(:,:,cropidx1:cropidx2);
    Mz0B_MAP = Mz0B_MAP(:,:,cropidx1:cropidx2);
    w = weights(:,:,cropidx1:cropidx2);
    
    % resample/downsample maps to desired size in x/y
    kTRANS = imresize3(kTRANS, matSize);
    k_1_2_MAP = imresize3(k_1_2_MAP, matSize);
    k_1_3_MAP = imresize3(k_1_3_MAP, matSize);
    Mz0P_MAP = imresize3(Mz0P_MAP, matSize);
    Mz0L_MAP = imresize3(Mz0L_MAP, matSize);
    Mz0B_MAP = imresize3(Mz0B_MAP, matSize);
    w = imresize3(w, matSize);
   
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
    
    kTRANS = flip(kTRANS,3);
    kMaps = flip(cat(4,k_1_2_MAP,k_1_3_MAP),3); %flip 3rd dimension to match in vivo convention
    Mz0Maps = flip(cat(4,Mz0P_MAP,Mz0L_MAP, Mz0B_MAP),3);
    w = flip(w,3);

    if ~isempty(simParams) % output simulated metabolite dynamic images
        % simulate signals
        % store simulation parameters a in a struct,
        % eventually this should be a custom class

        nMets = size(kineticRates,1) + 1;

        if isempty(inputFunction)
            inputFunction = zeros([1 simParams.Nt]);
        end
    
        metImages = zeros(cat(2,matSize,[nMets, simParams.Nt]));
        for Ix = 1:matSize(1)
            for Iy = 1:matSize(2)
                for Iz = 1:matSize(3)
                    Mz0_vx = [Mz0Maps(Ix, Iy, Iz, 1) Mz0Maps(Ix, Iy, Iz, 2) Mz0Maps(Ix, Iy, Iz, 3)];
                    [Mxy, ~] = simulate_Nsite_model(Mz0_vx, simParams.R1, [kMaps(Ix,Iy,Iz,1) 0; kMaps(Ix,Iy,Iz,2) 0], simParams.flips, simParams.TR, inputFunction*kTRANS(Ix,Iy,Iz) );
                    noise_S = randn([nMets simParams.Nt])* simParams.std_noise; % add noise %TO DO: define noise as percentage of input?
                    metImages(Ix,Iy,Iz,:,:) = Mxy + noise_S;
                end
            end
        end
        metImages = metImages .* repmat(w, [1 1 1 size(metImages,4) size(metImages,5)]);
    
    else
        metImages = 0;
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