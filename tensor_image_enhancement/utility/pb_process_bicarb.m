function [met_specAUC_bicarb,varargout] = pb_process_bicarb(filepath,params,Data_5D_pyrlac,varargin)

    do_TSE = 0;
    if nargin > 3,
        do_TSE = 1;
        sizeLR = varargin{1};
    end
   
    peaksearch_range = params.spec_idx_met{1};
    spec_idx_bicarb_new = {1:length(params.spec_idx_met{3})};
    load(filepath);
    % ---- pre-processing ----
    csi_phased = noise_prewhitening(csi_phased, params.dmtx);
    Data_5D_bicarb = squeeze(permute(csi_phased,[4 1 3 2 5 6]));
    Data_5D_bicarb = pb_correct_chemshift_brain(Data_5D_bicarb,peaksearch_range,params.circ_shift_dist);
    Data_5D_bicarb = Data_5D_bicarb .* ...
        repmat(permute(params.svdWeights,[4 2 3 1 5]),[size(Data_5D_bicarb,1) 1 1 1 size(Data_5D_bicarb,5)]);
    Data_5D_bicarb = squeeze(sum(Data_5D_bicarb,4));
    Data_5D_bicarb = pb_phase_baseline(Data_5D_bicarb,params.spec_idx_met);
    % ---- spectral splicing ----
    Data_5D_splice = pb_splicing(Data_5D_pyrlac,Data_5D_bicarb,params.spec_idx_met);
    % ---- TSE ----
    if do_TSE,
        % Tensor thresholding
        Data_5D_splice_LR = td_regularization(Data_5D_splice,sizeLR,0);
    end

    met_specAUC_before = spec2img(Data_5D_splice,spec_idx_bicarb_new);
    met_specAUC_after = spec2img(Data_5D_splice_LR,spec_idx_bicarb_new);
    
    met_specAUC_bicarb{1} = met_specAUC_before{1};
    met_specAUC_bicarb{2} = met_specAUC_after{1};
    
    varargout{1} = Data_5D_splice_LR;
    
end