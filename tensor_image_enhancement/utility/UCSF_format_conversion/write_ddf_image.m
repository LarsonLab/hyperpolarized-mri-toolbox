function write_ddf_image(filename, data, compressflag, displayflag)
%
%   WRITE_DDF_IMAGE Writes a DDF file and associated complex file
%
%   WRITE_DDF_IMAGE(ROOTNAME, DATA) writes the data contained
%   in DATA into ROOTNAME.DDF and ROOTNAME.CMPLX.  DATA should be
%   the result of a READ_DDF_IMAGE operation (though, of
%   course, you are allowed to modify the various parameters
%   contained in the DATA structures.)  The DATA can be v5 or v6.
%
%   WRITE_DDF_IMAGE(ROOTNAME, DATA, COMPRESSFLAG) will gzip the
%   data if COMPRESSFLAG == 1.
%
%   WRITE_DDF_IMAGE(ROOTNAME, DATA, COMPRESSFLAG, DISPLAYFLAG) will
%   do the same thing, but if DISPLAYFLAG == 0, it will not write
%   any information to the screen.
%
%   The recommended usage is to first read in some data set using
%   READ_DDF_IMAGE, perform your processing in Matlab, manually
%   modify any parts of the DDF structure that you need to change
%   (especially important for Matlab are things like number of
%   points in each dimension), then write using this program.  Even
%   if you want to create a complex data set completely anew in
%   Matlab (for example, a simulated spectral dataset), it is MUCH
%   easier to find some DDF that is more or less like the data you
%   are creating, and modify it to look like what you want.  It is
%   more trouble than it is worth to try to create a DDF manually.
%
%   See also READDDF_FILE, WRITEDDF_FILE, WRITE_DDF_IMAGE.
%
%     Michael C. Lee, Ph.D.
%     Department of Radiology
%     University of California, San Francisco
%
%
%   Copyright (c) 2009 Regents of the University of California.
%   All rights reserved.  This software provided pursuant to a
%   license agreement containing restrictions on its disclosure,
%   duplication and use.  This notice must be embedded in or
%   attached to all copies, including partial copies, of the
%   software or any revisions or derivations thereof.
%
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/write_ddf_image.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%



if (nargin < 2)
  error(['Not enough input arguments: must specify at least rootname '...
         'and a data structure']);
end

filename = char(filename);

if (nargin < 3) 
  compressflag = 0;
end

if (nargin < 4) 
  displayflag = 1;
end

if (displayflag == 1)
  disp(sprintf('\n  Writing DDF file:    %s.ddf',filename));
end

data.ddf.filename = filename;
data.ddf.rootname = strip_path(filename);

datafilename = strcat(filename,'.cmplx');

if (isfield(data,'img'))
  data_size = size(data.img);
else
  data_size = size(data.real);
end

if (ndims(data.img) == 4)
  data_size = data_size(2:4);
elseif (ndims(data.img) == 3)
  data_size = [data_size(2:3) 1];
elseif (ndims(data.img) == 2)
  data_size = [data_size(2) 1 1];
elseif (ndims(data.img) == 1)
  data_size = [1 1 1];
end

if (data_size(:) ~= data.ddf.npix(:));
  error(sprintf('\n-- DDF size: %d %d %d   Data size: %d %d %d! --\n\n', ...
	       data.ddf.npix, size(data.img)));
end

writeddf_file(filename, data.ddf);

datafile = fopen(datafilename,'wb','b');

if (displayflag == 1)
  disp(sprintf('  Writing DDF image:   %s',datafilename));
end

if (isfield(data,'img'))

  if (isstruct(data.img(1)))

    % Mike Style
    
    for (k = 1:data.ddf.npix(3))
      for (j = 1:data.ddf.npix(2))
        for (i = 1:data.ddf.npix(1)) 
          
          if (length(data.img(i,j,k).real) == 0)
            real_part = zeros(data.ddf.specpoints,1);
          else
            real_part = data.img(i,j,k).real;
          end
          
          if (length(data.img(i,j,k).imaginary) == 0)
            imaginary_part = zeros(data.ddf.specpoints,1);
          else
            imaginary_part = data.img(i,j,k).imaginary;
          end      
          
          this_vector(1:2:data.ddf.specpoints*2) = real_part;
          this_vector(2:2:data.ddf.specpoints*2) = imaginary_part;
          
          fwrite(datafile, this_vector, 'float');
          
        end    
      end
    end
  else
    
    % Esin Style
    
    data_vector = zeros(length(data.img(:))*2,1);
    data.img = data.img(:);
    data_vector(1:2:end) = real(data.img);
    data_vector(2:2:end) = imag(data.img);  
    fwrite(datafile, data_vector, 'float');
  end

else
  
  % Janine Style

  data.img = complex(data.real, data.imag);
  data_vector = zeros(length(data.img(:))*2,1);
  data.img = data.img(:);
  data_vector(1:2:end) = real(data.img);
  data_vector(2:2:end) = imag(data.img);  
  fwrite(datafile, data_vector, 'float');  
  
end

  
if (compressflag == 1)
  unix(sprintf('gzip -f -q %s', datafilename));
end

fclose(datafile);
