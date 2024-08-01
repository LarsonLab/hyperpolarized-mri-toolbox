function [kTRANS, kMaps, Mz0Maps, metImages] = cardiac_metabolic_phantom(kineticRates, ktransScales, Mz0_constants, matSize, simParams, heart_idx)
% CARDIAC_METABOLIC_PHANTOM generates standardized 3-dimensional perfusion
%   and metabolism maps for simulated experiments. Supports 3 chemical pool
%   kinetic rate mapping.
%
%   Toolboxes required: Image Processing
%
%   Parameters:
%       kineticRates    = the kinetic rates to simulate, [(# of chemical
%                         pools) - 1, 4 (tissue types: LV, RV, LV Myocardium, RV Myocardium)]
%       ktransScales    = volume transfer constants for different tissue
%                         compartments [LV, RV, LV Myocardium, Rv Myocardium]
%       Mz0_constants   = Initial magnetization constants per metabolite per compartment,
%                         size [Nmets 3] in order LV, RV, LV Myocardium, RV Myocardium
%       matSize         = 1x3 vector for desired matrix size of each
%                         dimension 
%                         [met dim(nx,ny,nz)], default = [32 32 5]
%       simParams       = parameters used for the kinetic simulations: 
%                         Tarrival, Tbolus, TR, Nt, R1, flips; if empty won't
%                         generate metImages, default = empty struct
%       inputFunction   = [1 Nt] vector as input into each voxel
%
%   Optional Parameters:
%       heart_idx       = heart mask used for (1-20). default = random mask
%
%   Outputs:
%       kTRANS      = generated perfusion map
%       kMaps       = generated rate maps for 1->2 and 1->3
%       Mz0Maps     = generated Mz0 maps
%       metImages   = simulated metabolite dynamic images
%
% Copyright, 2024

    % parse input arguments
    arguments % !! some may not be robust enough after changes !!
        kineticRates (:,4) double {mustBeNumeric} = [0.1, 0.2, 0.3; 0, 0, 0]
        ktransScales (:,4) double {mustBeNumeric} = [1, 0.3, 0.3]
        Mz0_constants (:,4) double {mustBeNumeric} = [0, 0, 0]
        matSize (:,3) double {mustBeInteger} = [32, 32, 5]
        simParams struct = struct([])
        heart_idx {mustBeInteger, mustBePositive, mustBeNonzero} = []
    end

    nMets = size(kineticRates,1) + 1;

    % load masks
    if isempty(heart_idx)
        heart_idx = randi(20);
    end
    current_path = pwd;
    mask_path = fullfile(current_path,'util',num2str(heart_idx), 'cardiac_masks.mat');   
    load(mask_path,'im_mask'); 

    nTissues = 4;
    if size(im_mask,4) ~= nTissues
        error('Unexpected number of tissues present in the imported masks, ask for help, idk.');
    end
    
    % parameters for output map generation
    permuted_mask = double(permute(im_mask,[4 1 2 3]));
    sumWeights = sum(im_mask,4);

    kTRANS = create_map(permuted_mask, ktransScales, sumWeights);

    % generate Tarrival maps
    Tarrival_MAP = create_map(permuted_mask, simParams.Tarrival, sumWeights);
    Tarrival_MAP(sumWeights == 0) = (simParams.Tarrival(3)+simParams.Tarrival(4)) / 2;
    Tarrival_MAP = squeeze(Tarrival_MAP);

    % generate Mz0
    % this is kinda awkward because you need to regenerate the input
    % functions later
    Mz0 = zeros(nMets, nTissues);
    for tissue = 1:nTissues
        input_function = realistic_input_function(simParams.Nt, simParams.TR, simParams.Tarrival(tissue), simParams.Tbolus);
        Mz0(:,tissue) = input_function(1) .* Mz0_constants(:,tissue);
    end
        
    % generate the kinetic rate maps
    k_1_2_MAP = create_map(permuted_mask, kineticRates(1,:), sumWeights);
    k_1_3_MAP = create_map(permuted_mask, kineticRates(2,:), sumWeights);

    % generate Mz0 maps
    Mz0P_MAP = create_map(permuted_mask, Mz0(1,:), sumWeights);
    Mz0L_MAP = create_map(permuted_mask, Mz0(2,:), sumWeights);
    Mz0B_MAP = create_map(permuted_mask, Mz0(3,:), sumWeights);
    
    small_matSize = [matSize(1)/2, matSize(2)/2, matSize(3)];

    % resample/downsample maps to desired matrix size
    kTRANS = imresize3(kTRANS, small_matSize);
    k_1_2_MAP = imresize3(k_1_2_MAP, small_matSize);
    k_1_3_MAP = imresize3(k_1_3_MAP, small_matSize);
    Mz0P_MAP = imresize3(Mz0P_MAP, small_matSize);
    Mz0L_MAP = imresize3(Mz0L_MAP, small_matSize);
    Mz0B_MAP = imresize3(Mz0B_MAP, small_matSize);
    Tarrival_MAP = imresize3(Tarrival_MAP, small_matSize, "linear");
    
    % upsample maps for blur
    kTRANS = imresize3(kTRANS, matSize);
    k_1_2_MAP = imresize3(k_1_2_MAP, matSize);
    k_1_3_MAP = imresize3(k_1_3_MAP, matSize);
    Mz0P_MAP = imresize3(Mz0P_MAP, matSize);
    Mz0L_MAP = imresize3(Mz0L_MAP, matSize);
    Mz0B_MAP = imresize3(Mz0B_MAP, matSize);
    Tarrival_MAP = imresize3(Tarrival_MAP, matSize, "linear");
    
    % account for downsampling changing kinetic rates
    % TODO: not a perfect solution, areas that should have lower rates go higher
    k_1_2_resize_error = max(kineticRates(1,:),[],'all') / max(k_1_2_MAP,[],'all');
    k_1_3_resize_error = max(kineticRates(2,:),[],'all') / max(k_1_3_MAP,[],'all');
    kTRANS_resize_error = max(ktransScales,[],'all') / max(kTRANS,[],'all');
    k_1_2_MAP = k_1_2_MAP .* k_1_2_resize_error;
    k_1_3_MAP = k_1_3_MAP .* k_1_3_resize_error;
    kTRANS = kTRANS .* kTRANS_resize_error;


    kMaps = cat(4,k_1_2_MAP,k_1_3_MAP); 
    Mz0Maps = cat(4,Mz0P_MAP,Mz0L_MAP, Mz0B_MAP);

    if ~isempty(simParams) % output simulated metabolite dynamic images
        % simulate signals
        % store simulation parameters a in a struct,
        
        % generate voxelwise dynamics
        metImages = zeros(cat(2,matSize,[nMets, simParams.Nt]));
        for Ix = 1:matSize(1)
            for Iy = 1:matSize(2)
                for Iz = 1:matSize(3)
                    % generate input_function
                    Tarrival_vx = Tarrival_MAP(Ix, Iy, Iz);
                    input_function_vx = realistic_input_function(simParams.Nt, simParams.TR, Tarrival_vx, simParams.Tbolus);
                    Mz0_vx = [Mz0Maps(Ix, Iy, Iz, 1) Mz0Maps(Ix, Iy, Iz, 2) Mz0Maps(Ix, Iy, Iz, 3)];
                    [metImages(Ix,Iy,Iz,:,:), ~] = simulate_Nsite_model(Mz0_vx, simParams.R1, [kMaps(Ix,Iy,Iz,1) 0; kMaps(Ix,Iy,Iz,2) 0], simParams.flips, simParams.TR, input_function_vx*kTRANS(Ix,Iy,Iz) );
                end
            end
        end
    end     
    
    % add rician noise
    metImages_w_noise = zeros(size(metImages));
    for Imet=1:nMets
        ImetImage = squeeze(metImages(:,:,:,Imet,:));
        std_noise = max(sum(ImetImage,4),[],'all') ./ (simParams.SNR(Imet) * sqrt(simParams.Nt));
        noise_R = randn(cat(2,matSize,[simParams.Nt]))* std_noise;
        noise_I = randn(cat(2,matSize,[simParams.Nt]))* std_noise;
        metImages_w_noise(:,:,:,Imet,:) = sqrt((ImetImage+ noise_R).^2 + noise_I.^2);
    end
    
    metImages = metImages_w_noise;

end


function [map] = create_map(mask, rates, sumWeights)
    % apply rates over masks
    map_wSum = pagemtimes(rates,mask);
    map = squeeze(map_wSum)./sumWeights;
    map(isnan(map)) = 0; 
end

