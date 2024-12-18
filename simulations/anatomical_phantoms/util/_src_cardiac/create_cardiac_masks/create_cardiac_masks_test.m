% testing script for create_cardiac_masks
clear; close all;
% PCA modes are from the Cardiac Atlas Project and was derived from the UK
% Biobank Study

% the HDF5 files are too large for GitHub, but can be found here:
% https://www.cardiacatlas.org/biventricular-modes/
% after download, place the .h5 file in the same directory as this script

% read point cloud data
sourceFile = './UKBRVLV_All.h5';
PCA_mode = 20; 
pc = h5read(sourceFile, '/COEFF');
ev = h5read(sourceFile, '/LATENT');
mu = h5read(sourceFile, '/MU');
S = mu + (1.5 .* sqrt(ev(1)) .* pc(:,PCA_mode)');
N = length(S);
ed = reshape(S(1:N/2), 3, [])';
es = reshape(S((N/2+1):end), 3, [])';

vertices = struct;
% define areas of the heart
% these was separated by hand
lv = ed(1:1428,:);
rv = ed(1429:3056,:);
mc = ed(3057:end,:);

figure('Color', 'w');
plot3(lv(:,1), lv(:,2), lv(:,3), 'b.');
hold on;
plot3(rv(:,1), rv(:,2), rv(:,3), 'r.');
plot3(mc(:,1), mc(:,2), mc(:,3), 'g.');
axis vis3d
axis equal

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