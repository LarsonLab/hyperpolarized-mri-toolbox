function pb_gen_dummy_mask(Data_5D,ref_mask_path)

idf_temp = read_idf_image(ref_mask_path);
idf_temp.img = ones(size(Data_5D,2),size(Data_5D,3));
idf_temp.idf.npix = [size(Data_5D,2) size(Data_5D,3) 1];
write_idf_image('dummy_mrs_mask',idf_temp);

end