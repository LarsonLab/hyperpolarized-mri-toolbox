# `unified_phantom`

A simple usage example can be found in `concise_unified_phantom_test.m`, and a more verbose example can be found in `verbose_unified_phantom_test.m`

`tissue_structure` - class that defines a structural model
- `obj.create_k_trans_map` - method that creates a kTRANS map

`metabolite` - class that defines properties of a metabolite in a tissue

`pk_params` - class that defines pharmacokinetic model parameters

`pk_model` - static class with methods to run a pharmacokinetic model
- `pk_model.run_pk_model` - wrapper that runs pharmacokinetic model with parameters defined as arrays
- `pk_model.generate_met_dynamics` - generates metabolite dynamics from pharmacokinetic parameters
- `pk_model.generate_all_met_dynamics` - generates array of metabolite dynamics from parameters defined as arrays
- `pk_model.generate_met_images` - generates metabolite dynamic images from met dynamics and a structural model
- `pk_model.apply_k_trans` - creates metabolite dynamic images from a set of low kTRANS images, a set of high kTRANS images, and a kTRANS map

`mri_system` - static class with methods to augment metabolite dynamic images
- `mri_system.run_mri_system` - wrapper that runs mri system model
- `mri_system.make_met_images_multires` - converts single-resolution metabolite images to multiresolution
- `mri_system.add_rician_noise` - adds Rician noise to multiresolution images
- `mri_system.upsample_to_output_size` - upsamples multiresolution images to desired output size
- `mri_system.augment` - applies augmentations to metabolite dynamic images
- `mri_system.apply_coil_lim` - adds coil limits to metabolite dynamice images
