function ddf = set_ddf_dyn_dim(ddf, name, ndyn)
%pcao 15.04.26 add if ndyn >1, determine if it is dynamic data 
if ndyn > 1
    ddf.filename = [name '.ddf'];
    ddf.root_name = name;
    ddf.numdim = ddf.numdim + 1;  %add one dim
    ddf.dimension_type{ddf.numdim} = 'unknown';    %dynamic direction
%     if ddf.npix
    ddf.npix(ddf.numdim-1) = ndyn;        %set the num of steps in this dim
    ddf.pixel_spacing = [ddf.pixel_spacing 0]; %set the pixel sapcing zero for this dim
end
end
