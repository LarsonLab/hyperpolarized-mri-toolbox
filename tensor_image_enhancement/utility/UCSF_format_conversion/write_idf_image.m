function [imagefile] = write_idf_image(filename, imagedata, dynam_pathname, compressflag, displayflag)
%
%   WRITE_DDF_IMAGE Writes a DDF file and associated complex file
%
%   WRITE_IDF_IMAGE(ROOTNAME, DATA) writes the data contained
%   in DATA into ROOTNAME.IDF and ROOTNAME.{EXT} where {EXT}
%   is .byt, .int2, or .real, according to the FILETYPE given in 
%   the IDF structure.  DATA is usually the result of a 
%   READ_IDF_IMAGE operation (though, of course, you are allowed 
%   to modify the various parameters contained in the DATA 
%   structures.)  The DATA is always written as v5.
%
%   WRITE_IDF_IMAGE(ROOTNAME, DATA, PATH) will check if the
%   path given in PATHNAME is writeable.  PATH does not 
%   automatically cause your file to be written to this path, it
%   merely returns an error if it finds that it can't write here.
%
%   WRITE_IDF_IMAGE(ROOTNAME, DATA, PATH, COMPRESSFLAG) will
%   gzip the output image if COMPRESSFLAG is 1.
%
%   WRITE_IDF_IMAGE(ROOTNAME, DATA, PATH, COMPRESSFLAG, DISPLAY)
%   If DISPLAYFLAG == 0, it will not write any information to the 
%   screen (ie. 'quiet' mode)
%
%   The recommended usage is to first read in some data set using
%   READ_IDF_IMAGE, perform your processing in Matlab, manually
%   modify any parts of the IDF structure that you need to change
%   (especially important for Matlab are things like number of
%   points in each dimension), then write using this program.  Even
%   if you want to create an image completely anew in Matlab
%   (for example, a simulated spectral dataset), it is MUCH
%   easier to find some IDF that is more or less like the data you
%   are creating, and modify it to look like what you want.  It is
%   more trouble than it is worth to try to create a IDF manually.
%
%   See also READIDF_FILE, WRITEIDF_FILE, WRITE_IDF_IMAGE.
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
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/write_idf_image.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%



filename = char(filename);

if (nargin < 4) 
  compressflag = 0;
end

if (nargin < 5) 
  displayflag = 1;
end

%if (displayflag == 1)
%  disp(sprintf('\n  Writing IDF file:    %s.idf',filename));
%end

switch imagedata.idf.filetype
 case 2
  extension = '.byt';
  datatype  = 'uchar';
 case 3
  extension = '.int2';
  datatype  = 'int16';
 case 7
  extension = '.real';
  datatype  = 'float32';
 case 17
  extension = '.cmplx';
  datatype  = 'float32';
 otherwise
  error('Unrecognized file type "%d"', ...
        output.idf.filetype);
end

imagedata.idf.filename = filename;

imagefilename = strcat(filename,extension);

if (ndims(imagedata.img) == 3)
  if (any(size(imagedata.img) ~= imagedata.idf.npix))
    error('size mismatch: image is %d %d %d, but IDF claims %d %d %d', ...
          size(imagedata.img,1), size(imagedata.img,2), ...
          size(imagedata.img,3), imagedata.idf.npix(1), ...
          imagedata.idf.npix(2), imagedata.idf.npix(3));
  end
else
  if (any(size(imagedata.img) ~= imagedata.idf.npix(1:2)));
    error('size mismatch: image is %d %d %d, but IDF claims %d %d %d', ...
          size(imagedata.img,1), size(imagedata.img,2), ...
          size(imagedata.img,3), imagedata.idf.npix(1), ...
          imagedata.idf.npix(2), imagedata.idf.npix(3));
  end
end

imagedata.idf.minimum = min([0 min(min(min(imagedata.img)))]);
imagedata.idf.maximum = max(max(max(imagedata.img)));

imagefile= fopen(imagefilename,'wb','b');

if imagefile == -1
  if (nargin > 2)
    error('Do not have access to write to the directory "%s"', dynam_pathname);
  else
    error('Unable to open file "%s" for writing', imagefilename);  
  end
else
   
   if (displayflag == 1)
     disp(sprintf('  Writing IDF image:   %s',imagefilename));
   end

   if (imagedata.idf.filetype ~= 17)   
     fwrite(imagefile, imagedata.img, datatype);
   else     
     temp_data(1:2:2*prod(imagedata.idf.npix)) = real(imagedata.img);
     temp_data(2:2:2*prod(imagedata.idf.npix)) = imag(imagedata.img);
     fwrite(imagefile, temp_data, datatype);
   end
   
   if (compressflag == 1)
     unix(sprintf('gzip -f -q %s', imagefilename));
   end
   
   fclose(imagefile);

   writeidf_file(filename, imagedata.idf)
   
   if (displayflag == 1)
     disp(sprintf('  Writing IDF file:    %s.idf',filename));
   end 

end

