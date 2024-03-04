% Quick testing script for the brain web basic metabolic phantom
kineticRates = [0, 0.04, 0.02; 0, 0, 0]; % no kPB, effectively, in order [vasc, GM, WM]
ktransScales = [1, 0.3 0.3];
isFuzzy = true;

[k_trans, k_maps] = brainweb_metabolic_phantom(kineticRates,ktransScales,isFuzzy);

%% check a few central slices

%slices = 100:10:200;
slices = 1:8;

% visualize kTRANS
figure,
imagescn(k_trans(:,:,slices),[0 max(k_trans(:,:,slices),[],'all')], [1 numel(slices)]); colormap fire;

% visualize kinetic rates
figure,
imagescn(k_maps(:,:,slices,1),[0 max(k_maps(:,:,slices,1),[],'all')], [1 numel(slices)]); colormap fire;

figure, % picking an arbitrarily small number for the lower scale bound since k13 == 0
imagescn(k_maps(:,:,slices,2),[-1e-12 1e-12], [1 numel(slices)]); colormap fire;

