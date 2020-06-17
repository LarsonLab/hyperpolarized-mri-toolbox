function [dataOut,varargout] = pb_correct_chemshift_brain(dataIn,peaksearch_range,varargin)
% correct B0 spectral shifts
% peaksearch_range: index range of pyruvate peak in ref voxel
% dataIn(f,x1,x2,x3...,dyn)
calc_B0_flag = 1;
if nargin > 2
    calc_B0_flag = 0;
    circ_shift_dist = varargin{1}(:);
    fprintf('Correct chemshift based on table ... \n');
end

dataIn_AUC = sum(abs(dataIn),ndims(dataIn));
size_dataIn = size(dataIn);
num_voxel = prod(size_dataIn(2:end-1));
if calc_B0_flag
    dataIn_AUC_flat = reshape(dataIn_AUC,size_dataIn(1),[]);
    load('model_spectra.mat');
    spec_ref = model_spectra;
    [~,empirical_center_temp] = max(spec_ref(peaksearch_range));
    empirical_center = empirical_center_temp+peaksearch_range(1)-1;
    circ_shift_dist = zeros(1,num_voxel);

    notfound_cnt = 0;
    for i_voxel = 1:num_voxel,
        [notfound closest_Idx] = peak_chemshift_corr(spec_ref,...
            squeeze(dataIn_AUC_flat(:,i_voxel))',...
            empirical_center,64,30);
        circ_shift_dist(i_voxel) = empirical_center-closest_Idx;
        notfound_cnt = notfound_cnt+notfound;
    end
    fprintf('cannot find shift in %01d voxels \n',notfound_cnt);
end

dataIn_temp = reshape(dataIn,size_dataIn(1),num_voxel,[]);
for i_voxel = 1:num_voxel,
for i_dyn = 1:size_dataIn(end),
    dataOut(:,i_voxel,i_dyn) = circshift(dataIn_temp(:,i_voxel,i_dyn),circ_shift_dist(i_voxel));
end
end
dataOut = reshape(dataOut,size_dataIn);
varargout{1} = reshape(circ_shift_dist,size_dataIn(2:end-1));

end