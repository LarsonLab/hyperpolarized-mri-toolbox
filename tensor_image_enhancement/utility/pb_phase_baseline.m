function Data_5D_Out = pb_phase_baseline(Data_5D_In,spec_idx_met,varargin)

    % ---- low order phase ----
    Data_5D_In = pb_phase_corr(Data_5D_In,spec_idx_met);
    
    % ---- Baseline Correction ----
    peakfile_name = '2d_epsi_c13';
    if nargin > 2
        peakfile_name = varargin{1};
    end
    % save
    dir_temp = 'EPSI_spec_temp';
    dummypath = 'utility/reference_spectra/dynout_phased_1_01';
    flip_dims = [2];
    pb_saveddf(Data_5D_In,dir_temp,dummypath,flip_dims);
    % baseline
    pb_gen_dummy_mask(Data_5D_In,'t10798_dynout_mrs_mask');
    i_dyn = 1;
    while(exist(sprintf('%s/dynaa_timepoint%02d.ddf',dir_temp,i_dyn),'file'))
        system(sprintf('c13_epsi_baseline_correction_v2.x %s/dynaa_timepoint%02d %s',...
            dir_temp,i_dyn,peakfile_name));
        i_dyn = i_dyn +1;
    end
    delete('dummy_mrs_mask.*')
    % reload corrected spectrum
    Data_5D_Out = pb_load_spectrum(sprintf('%s/dynaa_timepoint',dir_temp),flip_dims);
    % delete temp
    rmdir(dir_temp,'s');

end