# `unified_phantom`

A simple usage example can be found in `concise_unified_phantom_test.m`, and a more verbose example can be found in `verbose_unified_phantom_test.m`

`tissue_structure` - class that defines a structural model

`metabolite` - class that defines properties of a metabolite in a tissue

`pk_params` - class that defines pharmacokinetic model parameters

`pk_model` - static class with methods to run a pharmacokinetic model
- `pk_model.run_pk_model` - runs pharmacokinetic model with parameters defined as arrays
- `pk_model.generate_met_dynamics` - generates metabolite dynamics from pharmacokinetic parameters
- `pk_model.generate_met_images` - generates metabolite dynamic images from met dynamics and a structural model

`mri_system` - static class with methods to augment metabolite dynamic images
- `mri_system.run_mri_system` - runs mri system model with parameters defined as arrays
- `mri_system.make_met_images_multires` - converts single-resolution metabolite images to multiresolution
- `mri_system.add_rician_noise` - adds Rician noise to multiresolution images
