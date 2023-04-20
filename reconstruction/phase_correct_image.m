function image_corrected = phase_correct_image(image_original)
% image_corrected = phase_correct_image(image_original)
%
% perform voxel-wise phase correction, aiming to put signal in the real
% channel
%
% image_original - final dimension is the dimension that phase correction
% is applied along, e.g. should be [x,y,t], [x,y,f], [x,y,z,t], [x,y,z,f]


size_image = size(image_original);

image_original = reshape(image_original, [prod(size_image(1:end-1)), size_image(end)]);
image_corrected = zeros(size(image_original));

for n = 1:size(image_original,1)
    image_corrected(n,:) = image_original(n,:)* exp(1i*find_phase_corr(image_original(n,:)));
end

image_corrected = reshape(image_corrected, size_image);