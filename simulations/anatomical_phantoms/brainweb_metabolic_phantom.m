function [kTRANS, kMaps, metImages] = brainweb_metabolic_phantom(kineticRates, ktransScales, isFuzzy, nAugment, augmentSeed)%, imagesOutFlag)
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
%       nAugment        = number of random augmentations, if not defined a
%                         single unaugmented phantom will be generated,
%                         default = 1
%       augmentSeed     = random seed for augmentations
%       imagesOutFlag   = flag indicating whether metabolite images will be
%                         simulated and returned as an output, default = 0
%
%   Outputs:
%       kTRANS      = generated perfusion map
%       kMaps       = generated rate maps for 1->2 and 1->3
%       metImages   = (UNDER CONSTRUCTION)simulated metabolite dynamic images
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
    nAugment double {mustBeInteger, mustBePositive, mustBeNonzero} = 1
    augmentSeed double {mustBeInteger, mustBePositive, mustBeNonzero} = 1
    %imagesOutFlag double {mustBeNumericOrLogical} = false
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
kMaps = cat(4,k_1_2_MAP,k_1_3_MAP);
kMaps(isnan(kMaps)) = 0;


% Add Augmentation support here eventually
% have fun Sule :)

metImages = 0;
% UNDER CONSTRUCTION: output of simulated images is not completed
% if imagesOutFlag % output simulated metabolite dynamic images
%     % simulate signals
%     % assume simulation parameters are stored in a cell array of structs,
%     % because I said so, eventually this should be a custom class
%     params = simParams{1}; % Mz0, Tarrival, Tbolus, TR, Nt, R1,
% 
%     input_function = realistic_input_function(params.Nt, params.TR, params.Tarrival, params.Tbolus);% normalize for a total magnetization input = 1
%     Mz0 = [0,0];
%     nMets = size(kineticRates,1) + 1;
% 
%     for i=1:numel(nTissues)
%         [Mxy, ~] = simulate_Nsite_model(Mz0, params.R1, [k_1_2(i) k_1_3(i)], params.flips, params.TR, input_function*kTRANS(k) );
%         noise_S = randn([2 params.Nt])* params.std_noise;
%         simData{i} = Mxy + noise_S;
%     end
% 
%     % create metabolite masks per tissue type
%     % vascular pool will reflect only the input function for now
%     metImages = zeros([size(vasc_mask) nMets]);
% 
% end

end