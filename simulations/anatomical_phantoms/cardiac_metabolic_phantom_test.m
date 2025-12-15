% quick testing script for cardiac_metabolic_phantom
clear; 
close all;

%% Arguments
kTRANS_scales = [1 1 -0.5 -0.5]; 
sample_size = [32 32 11; 
               16 16 11; 
               24 24 11];
output_size = [64 64 11];
kinetic_rates = [0.0075, 0.0045, 0.06, 0.02;
                0.0011, 0.0005, 0.0400, 0.01];

sim_params.Tarrivals = [7 0 10 14];
sim_params.Tbolus = 1;
sim_params.TR = 3.6;
sim_params.Nt = 30;
sim_params.R1 = [1/30 1/25 1/25];
sim_params.flips = repmat([20; 30; 30],[1 sim_params.Nt])*pi/180;
sim_params.SNR = [150 40 20]; 
sim_params.coil_lim = [0.4 1.2];

[input_functions, Mz0] = generate_input_functions(sim_params, 4);


function [input_functions, Mz0] = generate_input_functions(sim_params, n_tissues)
    input_functions = zeros(4, sim_params.Nt);
    for i = 1:n_tissues
        input_functions(i,:) = realistic_input_function(sim_params.Nt, sim_params.TR, sim_params.Tarrivals(i), sim_params.Tbolus);
    end

    Mz0 = [input_functions(1,1), input_functions(2,1), input_functions(3,1)*.5, input_functions(4,1)*.5;
           0, 0, input_functions(3,1)*.01, input_functions(4,1)*.01;
           0, 0, input_functions(3,1)*.005, input_functions(4,1)*.005;];
end




% Phantom
[kinetic_maps, kTRANS, Mz0_maps, input_function_map, met_images] ...
= cardiac_metabolic_phantom(kinetic_rates, kTRANS_scales, Mz0, input_functions, sample_size, output_size, sim_params);

slices = 1:size(kTRANS,3);

%% kTRANS
kTRANS_sum = sum(kTRANS, 4);
figure("Name","kTRANS")
imagescn(kTRANS_sum(:,:,slices), ...
    [0 max(kTRANS_sum(:,:,slices), [], 'all')], ...
    [1 numel(slices)]); 
colormap hot;

%% Kinetic maps
kinetic_sum = sum(kinetic_maps, 4);
figure("Name","kinetic 1->2")
imagescn(kinetic_sum(:,:,slices,1), ...
    [0 max(kinetic_sum(:,:,slices,1), [], 'all')], ...
    [1 numel(slices)]); 
colormap hot;

figure("Name","kinetic 1->3")
imagescn(kinetic_sum(:,:,slices,2), ...
    [0 max(kinetic_sum(:,:,slices,2), [], 'all')], ...
    [1 numel(slices)]); 
colormap hot;

%% Mz0 maps
Mz0_sum = sum(Mz0_maps, 4);
figure("Name","Mz0 Pyruvate")
if (max(Mz0_sum(:,:,slices,1), [], 'all')) > 0
    imagescn(Mz0_sum(:,:,slices,1), ...
        [0 max(Mz0_sum(:,:,slices,1), [], 'all')], ...
        [1 numel(slices)]); 
    colormap hot;
end

if (max(Mz0_sum(:,:,slices,2), [], 'all')) > 0
    figure("Name","Mz0 Lactate")
    imagescn(Mz0_sum(:,:,slices,2), ...
        [0 max(Mz0_sum(:,:,slices,2), [], 'all')], ...
        [1 numel(slices)]); 
    colormap hot;
end

if (max(Mz0_sum(:,:,slices,3), [], 'all')) > 0
    figure("Name","Mz0 Bicarb")
    imagescn(Mz0_sum(:,:,slices,3), ...
        [0 max(Mz0_sum(:,:,slices,3), [], 'all')], ...
        [1 numel(slices)]); 
    colormap hot;
end

%% met images
time_pts = 1:20;
figure("Name", "Pyruvate")
imagescn(met_images{1}(:,:,5,time_pts), ...
    [0 max(met_images{1}(:,:,5,time_pts), [], 'all') / 3] ...
    );
colormap hot;

figure("Name", "Lactate")
imagescn(met_images{2}(:,:,5,time_pts), ...
    [0 max(met_images{1}(:,:,5,time_pts), [], 'all') / 3] ...
    );
colormap hot;

figure("Name", "Bicarb")
imagescn(met_images{3}(:,:,5,time_pts), ...
    [0 max(met_images{1}(:,:,5,time_pts), [], 'all') / 3] ...
    );
colormap hot;
