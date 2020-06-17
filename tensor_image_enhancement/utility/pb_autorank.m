function [sizeLR_auto, varargout] = pb_autorank(Data_5D_orig,spec_idx_met,varargin)

    BVratio = 1; % bia-variance weighting in autorank, default = 1
    if nargin > 2
        BVratio = varargin{1}; % use when downsampled back to raw resolution
    end
    if nargin > 3
        tradeOff_flag = true;
    end
    % set up for auto rank calculations
    size_Data_5D_orig = size(Data_5D_orig);
    rankSweepRange = {5:20,1:size_Data_5D_orig(2),1:size_Data_5D_orig(3),1:size_Data_5D_orig(4)};
    imgIn = Data_5D_orig;
    [U,S,sv] = mlsvd(imgIn);
    ndim_imgIn = ndims(imgIn);
    imgIn_temp.U = U;
    imgIn_temp.S = S;
    imgIn_temp.sv = sv;
    imgIn_temp.ndim_imgIn = ndim_imgIn;

    spec_idx_met_pyrlac = spec_idx_met(1:2);
    
    TdVariance = NaN(size_Data_5D_orig);
    TdBias = NaN(size_Data_5D_orig);
    
    poolobj = parpool;
    parfor i_spec = rankSweepRange{1},
        fprintf('starting spec %d \n',i_spec);
        temp_TdVariance = NaN(size_Data_5D_orig(2:4));
        temp_TdBias = NaN(size_Data_5D_orig(2:4));
        for i_x = rankSweepRange{2},
        for i_y = rankSweepRange{3},
        for i_dyn = rankSweepRange{4},
        fprintf('calculating rank = (%d,%d,%d,%d) \n',i_spec,i_x,i_y,i_dyn);
        % Tensor thresholding
        sizeLR = [i_spec, i_x, i_y, i_dyn];
        Data_5D_trunc = td_regularization_fast(imgIn_temp,sizeLR,0);

        % resonance peak AUC estimates
        met_specAUC_before = spec2img(flip(Data_5D_orig,3),spec_idx_met_pyrlac);
        met_specAUC_after = spec2img(flip(Data_5D_trunc,3),spec_idx_met_pyrlac);
        
        % estimate bias and variance
        temp_TdVariance(i_x,i_y,i_dyn) = td_EstimateVariance(Data_5D_trunc);
        temp_TdBias(i_x,i_y,i_dyn) = td_EstimateBias(Data_5D_orig, Data_5D_trunc,spec_idx_met_pyrlac);  

        end
        end
        end
        TdVariance(i_spec,:,:,:) = temp_TdVariance;
        TdBias(i_spec,:,:,:) = temp_TdBias;
    end
    delete(poolobj);

    [CF_Params,TdCostFxn] = td_CostFxnTrade(TdBias,TdVariance,BVratio);
    Td_Params.TdCostFxn = TdCostFxn;

    sizeLR_auto = CF_Params.SortedRank(CF_Params.IdxMinCost,:);
    Td_Params.TdBias = TdBias;
    Td_Params.TdVariance = TdVariance;

    % optional output parameters
    varargout{1} = CF_Params;
    varargout{2} = Td_Params;

end