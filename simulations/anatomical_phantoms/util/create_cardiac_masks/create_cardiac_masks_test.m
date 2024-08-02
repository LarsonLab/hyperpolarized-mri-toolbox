% testing script for create_cardiac_masks
clear; close all;
% PCA modes are from the Cardiac Atlas Project and was derived from the UK
% Biobank Study

% the HDF5 files are too large for GitHub, but can be found here:
% https://www.cardiacatlas.org/biventricular-modes/

% read point cloud data
PCA_mode = 20; 
pc = h5read('UKBRVLV_All.h5', '/COEFF');
ev = h5read('UKBRVLV_All.h5', '/LATENT');
mu = h5read('UKBRVLV_All.h5', '/MU');
S = mu + (1.5 .* sqrt(ev(1)) .* pc(:,PCA_mode)');
N = length(S);
ed = reshape(S(1:N/2), 3, [])';

vertices = struct;
% define areas of the heart
% these was separated by hand
lv = ed(1:1428,:);
rv = ed(1429:3056,:);
mc = ed(3057:end,:);

% put it all in a struct
vertices.mc = mc;
vertices.lv = lv;
vertices.rv = rv;

maskSize = [-100,100,201;
            -100,100,201;
            -50,50,11];

[im_mask, layered_masks] = create_cardiac_masks(vertices, 2, maskSize);

slices = 1:round(size(im_mask,3)/10):size(im_mask,3);

%% visualize im_mask

figure;
imagescn(im_mask(:,:,slices,:), [0, 1], [numel(slices), size(im_mask,4)]);

%% visualize layered_masks

figure;
imagescn(layered_masks(:,:,slices,:),[0 max(layered_masks(:,:,slices,:),[],'all')], [1, numel(slices)]);