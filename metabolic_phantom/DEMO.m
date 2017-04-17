%% generate 3D metabolic phantom at 16x16x16 resolution 

nx = 16; 
ny = 16;
nz = 16;
kTRANS_low  = 0.02;
kTRANS_high = 0.05; 
kPL_low     = 0.01;
kPL_high    = 0.03; 

[kTRANS, kPL] = metabolic_phantom(nx, ny, nz, kTRANS_low, kTRANS_high, kPL_low, kPL_high); 

% plot kTRANS and kPL maps 

z_slice = nz/2; % define slice to plot 

f = figure(1);
imagesc(kTRANS(:, :, z_slice))
xticks([])
yticks([])
title('k_{TRANS} map slice at z=0')
pbaspect([1 1 1])
colorbar()
saveas(f, 'kTRANS_16.png', 'png')

f = figure(2);
imagesc(kPL(:, :, z_slice))
xticks([])
yticks([])
title('k_{PL} map slice at z=0')
pbaspect([1 1 1])
colorbar()
saveas(f, 'kPL_16.png', 'png')


%% generate 3D metabolic phantom at 256x256x256 resolution 

nx = 256; 
ny = 256;
nz = 256;

[kTRANS, kPL] = metabolic_phantom(nx, ny, nz, kTRANS_low, kTRANS_high, kPL_low, kPL_high); 

% plot kTRANS and kPL maps 

z_slice = nz/2; % define slice to plot 

f = figure(3);
imagesc(kTRANS(:, :, z_slice))
xticks([])
yticks([])
title('k_{TRANS} map slice at z=0')
pbaspect([1 1 1])
colorbar()
saveas(f, 'kTRANS_256.png', 'png')

f = figure(4);
imagesc(kPL(:, :, z_slice))
xticks([])
yticks([])
title('k_{PL} map slice at z=0')
pbaspect([1 1 1])
colorbar()
saveas(f, 'kPL_256.png', 'png')
