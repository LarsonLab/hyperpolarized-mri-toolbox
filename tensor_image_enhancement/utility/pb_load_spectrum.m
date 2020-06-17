function dataOut = pb_load_spectrum(filepath,varargin)
    addpath('/netopt/share/lib/local/brain/matlab');
    flip_dims = [];
    if nargin > 1
        flip_dims = varargin{1};
    end
    num_files = 0;
    while(exist(sprintf('%s%02db_cor.ddf',filepath,num_files+1),'file'))
        ddf_temp = read_ddf_image(sprintf('%s%02db_cor',filepath,num_files+1));
        dataOut(:,:,:,num_files+1) = ddf_temp.img;
        num_files = num_files+1;
    end
    for i_dim = flip_dims,
        dataOut = flip(dataOut,i_dim);
    end
end