function [TdBias] = td_EstimateBias(Data_5D_before, Data_5D_after, spec_idx_met)
% v3: find top x across the first y timepoints
% ----------------------
% ImgBefore(x,y,dyn)
NumMaxDyn = 5; % use top x timepoints of SNR
NumMaxVoxels = 10; % use top x voxels of SNR

size_Data_5D = size(Data_5D_before);
met_specAUC_before = spec2img(flip(Data_5D_before,3),spec_idx_met);
met_specAUC_after = spec2img(flip(Data_5D_after,3),spec_idx_met);

dyn_idx = zeros(size_Data_5D(4),2);
TdBias = 0;

for i_met = 1:2,
    % now pick out top time points
    vol_dyn_temp = squeeze(sum(sum(met_specAUC_before{i_met},1),2));
    [~,dyn_idx(:,i_met)] = sort(vol_dyn_temp, 'descend');
    ImgBefore = met_specAUC_before{i_met};
    ImgAfter = met_specAUC_after{i_met};
    IdxMaxDyn = dyn_idx(1:NumMaxDyn,i_met)';
    % pick out largest voxels
    ImgBefore_AUC = sum(ImgBefore(:,:,IdxMaxDyn),3);
    ImgBefore_AUC = ImgBefore_AUC(:);
    [~,vx_idx] = sort(ImgBefore_AUC, 'descend');
    IdxMaxVoxel_timepoint = vx_idx(1:NumMaxVoxels);
    [IdxMaxVoxel_I,IdxMaxVoxel_J] = ind2sub(size_Data_5D(2:3),IdxMaxVoxel_timepoint);
    IdxMaxVoxel_K_temp = repmat(IdxMaxDyn,[NumMaxVoxels 1]);
    IdxMaxVoxel = sub2ind(size_Data_5D(2:4),repmat(IdxMaxVoxel_I,[NumMaxDyn 1]),...
        repmat(IdxMaxVoxel_J,[NumMaxDyn 1]),...
        IdxMaxVoxel_K_temp(:));
    % evaluate bias
    ValMaxVoxel_Before = ImgBefore(IdxMaxVoxel);
    ValMaxVoxel_After = ImgAfter(IdxMaxVoxel);
    TdBias_temp = mean(abs(ValMaxVoxel_After-ValMaxVoxel_Before));
    TdBias = TdBias+TdBias_temp;
end

end