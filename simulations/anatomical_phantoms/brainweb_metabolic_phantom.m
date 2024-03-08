function [kTRANS, kMaps, metImages] = brainweb_metabolic_phantom(kineticRates, ktransScales, isFuzzy, matSize, simParams, augmentParams, augmentSeed)
% BRAINWEB_METABOLIC_PHANTOM generates standardized 3-dimensional perfusion
%   and metabolism maps for simulated experiments. Supports 3 chemical pool
%   kinetic rate mapping.
%
%   Parameters:
%       kineticRates    = the kinetic rates to simulate, [(# of chemical
%                         pools) - 1, 3 (tissue types: vasc, gm, wm)]
%       ktransScales    = volume transfer constants for different tissue
%                         compartments [vasculature, gray matter, white matter]
%       isFuzzy         = logical flag to indicate if fuzzy tissue boundaries
%                         will be used, default = true
%       matSize         = 1x3 vector for desired matrix size of each dimension
%                         [nx ny nz], deafult = [16 16 8]
%       simParams       = parameters used for the kinetic simulations: Mz0, 
%                         Tarrival, Tbolus, TR, Nt, R1, flips; if empty won't
%                         generate metImages, default = empty struct
%       augmentParams   = random augmentation parameters XTranslation, 
%                         YTranslation etc, if not defined a single 
%                         unaugmented phantom will be generated,
%                         default = empty struct
%       augmentSeed     = random seed for augmentations, default no seed
%
%   Outputs:
%       kTRANS      = generated perfusion map
%       kMaps       = generated rate maps for 1->2 and 1->3
%       metImages   = simulated metabolite dynamic images
%
%   Author:
%       Jasmine Hu
%       Anna Bennett
%       Sule Sahin

    % parse input arguments
    arguments
        kineticRates (:,3) double {mustBeNumeric} = [0.1, 0.2, 0.3; 0, 0, 0]
        ktransScales (1,3) double {mustBeNumeric} = [1, 0.3, 0.3]
        isFuzzy double {mustBeNumericOrLogical} = true
        matSize (1,3) double {mustBeInteger} = [16, 16, 8]
        simParams struct = struct([])
        augmentParams struct = struct()
        augmentSeed double {mustBeInteger, mustBePositive, mustBeNonzero} = []
    end

    if augmentSeed
        rng(augmentSeed)
    else
        rng("shuffle")
    end
    
    % define default augmentation parameters and set augmentParams to them
    % if not defined
    defaultAugParams.XReflection = false; defaultAugParams.YReflection = false;
    defaultAugParams.Rotation = [0 0]; defaultAugParams.Scale = [1 1];
    defaultAugParams.XShear = [0 0]; defaultAugParams.YShear = [0 0];
    defaultAugParams.XTranslation = [0 0]; defaultAugParams.YTranslation = [0 0];
    augs = fields(defaultAugParams);
    for f=1:size(augs,1)
        if ~isfield(augmentParams,augs{f})
            augmentParams.(augs{f}) = defaultAugParams.(augs{f});
        end
    end

    % load base anatomical information
    resourcesDir = './resources';
    if isFuzzy
        brainwebFile = 'brainweb04_fuzzy_hires.mat';
    else
        brainwebFile = 'brainweb04_hi_res.mat';
    end
    baseMaskFile = dir(fullfile(resourcesDir,brainwebFile));
    if isempty(baseMaskFile)
        error('Error. \nBrainWeb mask file, %s, not found in resources directory.',brainwebFile)
    end
    
    load(fullfile(baseMaskFile.folder,brainwebFile),'im_mask');
    vasc_mask = squeeze(im_mask(:,:,1,:));
    maskSize = size(im_mask,1:3);
    gm_mask = squeeze(im_mask(:,:,2,:));
    wm_mask = squeeze(im_mask(:,:,3,:));
    
    k_1_2 = kineticRates(1,:);
    k_1_3 = kineticRates(2,:);
    
    nTissues = 3;
    if size(im_mask,4) ~= nTissues
        error('Unexpected number of tissues present in the imported masks, ask for help, idk.');
    end
    
    % parameters for output map generation
    permuted_mask = permute(im_mask,[4 1 2 3]);
    sumWeights = sum(im_mask,4);
    
    % generate the kTRANS masked volume
    kTRANS_wSum = pagemtimes(ktransScales,permuted_mask);
    kTRANS = squeeze(kTRANS_wSum)./sumWeights;
    kTRANS(isnan(kTRANS)) = 0;
    
    % generate the kinetic rate maps
    k_1_2_wSum = pagemtimes(k_1_2,permute(im_mask,[4 1 2 3]));
    k_1_2_MAP = squeeze(k_1_2_wSum)./sumWeights;
    k_1_3_wSum = pagemtimes(k_1_3,permuted_mask);
    k_1_3_MAP = squeeze(k_1_3_wSum)./sumWeights;
    k_1_2_MAP(isnan(k_1_2_MAP)) = 0;
    k_1_3_MAP(isnan(k_1_3_MAP)) = 0;
    
    % resample/downsample maps to desired size
    kTRANS = imresize3(kTRANS, matSize);
    k_1_2_MAP = imresize3(k_1_2_MAP, matSize);
    k_1_3_MAP = imresize3(k_1_3_MAP, matSize);
    
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
    
    kMaps = cat(4,k_1_2_MAP,k_1_3_MAP);
    
    if ~isempty(simParams) % output simulated metabolite dynamic images
        % simulate signals
        % store simulation parameters a in a struct,
        % eventually this should be a custom class
    
        input_function = realistic_input_function(simParams.Nt, simParams.TR, simParams.Tarrival, simParams.Tbolus);
        nMets = size(kineticRates,1) + 1;
    
        metImages = zeros(cat(2,matSize,[nMets, simParams.Nt]));
        for Ix = 1:matSize(1)
            for Iy = 1:matSize(2)
                for Iz = 1:matSize(3)
                    [Mxy, ~] = simulate_Nsite_model(simParams.Mz0, simParams.R1, [kMaps(Ix,Iy,Iz,1) 0; kMaps(Ix,Iy,Iz,2) 0], simParams.flips, simParams.TR, input_function*kTRANS(Ix,Iy,Iz) );
                    noise_S = randn([nMets simParams.Nt])* simParams.std_noise; % add noise %TO DO: define noise as percentage of input?
                    metImages(Ix,Iy,Iz,:,:) = Mxy + noise_S;
                end
            end
        end
    
    else
        metImages = 0;
    end

end