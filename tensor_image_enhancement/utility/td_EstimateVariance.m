function [TdVariance] = td_EstimateVariance(Data_5D)
% calculate across dyn points 01292019

NumLastDyn = 8;
% FID_noise_region = 2:12;
FID_noise_region = 8:16;

size_Data_5D = size(Data_5D);
for i_dim = 1:3,
    Data_5D = fft(fftshift(Data_5D,i_dim),[],i_dim);
end
noise_voxel = {[1 1],[1 size_Data_5D(3)],...
    [size_Data_5D(2) 1],[size_Data_5D(2) size_Data_5D(3)]};
noise_vector = [];
for i_voxel = 1:length(noise_voxel),
    temp = Data_5D(FID_noise_region,noise_voxel{i_voxel}(1),noise_voxel{i_voxel}(2),end-NumLastDyn:end);
    noise_vector = [noise_vector temp(:)'];
end
TdVariance = std(noise_vector);

end