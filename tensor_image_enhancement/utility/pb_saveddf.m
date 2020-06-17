function pb_saveddf(Data_5D,dir_output,dummypath,flip_dims)
    mkdir(dir_output);
    dummy_ddf = read_ddf_image(dummypath);
    dummy_ddf.ddf.specpoints = size(Data_5D,1);
    dummy_ddf.ddf.acq_n_data_points = size(Data_5D,1);
    dummy_ddf.ddf.npix = [size(Data_5D,2) size(Data_5D,3) 1];
    dummy_ddf.ddf.acq_n_points = dummy_ddf.ddf.npix;

    for i_dim = flip_dims,
        Data_5D = flip(Data_5D,i_dim);
    end
    for i_tp = 1:size(Data_5D,4),
        dummy_ddf.img = Data_5D(:,:,:,i_tp);
        write_ddf_image(sprintf('%s/dynaa_timepoint%02d',dir_output,i_tp),dummy_ddf);
    end
end