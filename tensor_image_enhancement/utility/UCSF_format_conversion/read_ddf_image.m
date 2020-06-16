function output = read_ddf_image(rootname, displayflag, complexformat)
%
%   READ_DDF_IMAGE Reads an DDF file and associated data file
%
%   DATA_STRUCTURE = READ_DDF_IMAGE(ROOTNAME) reads two files for this
%   dataset.  First, it reads 'ROOTNAME.ddf' which gives it the associated
%   header information.  This DDF can be either an MRSC v5 or v6 file.  
%   A data file is then read; these are always .cmplx files, consisting of
%   a set of 32-bit float pairs representing real and imaginary parts
%   of each value.  The results are returned as a structure with two fields.
%   The first is IMAGE_STRUCTURE.DDF where DDF is itself a structure with
%   fields corresponding to data entries found in the IDF file.  The second
%   is the complex data, which can come in a variety of user
%   specified formats (see below).
%
%   IMAGE_STRUCTURE = READ_DDF_IMAGE(ROOTNAME, DISPLAYFLAG) does the same
%   thing, but if DISPLAYFLAG is a 0, then it doesn't print out information
%   about what file it's opening.  This is useful in those cases when you
%   are trying to have a nice clean uncluttered screen.
%
%   IMAGE_STRUCTURE = READ_DDF_IMAGE(ROOTNAME, DISPLAYFLAG, FORMATFLAG) 
%   There are 3 integer options for the format of the complex data,
%   as described here:
%
%      0 (default)    a 4-D array of complex values, accessed as 
%                     IMAGE_STRUCTURE.IMG(FREQUENCY, X, Y, Z)
%                     ... which returns a complex value.
%
%      1              a 3-D array of structures.  Each element of
%                     the array is a structure with a REAL field
%                     and an IMAGINARY field.  Access these as
%                     IMAGE_STRUCTURE.IMG(X,Y,Z).REAL(FREQUENCY) and 
%                     ...IMAGINARY(FREQUENCY).
%                     (This is the old READ_DDF_IMAGE style)
%
%      2              Two 3-D arrays of doubles.  Access these as
%                     IMAGE_STRUCTURE.REAL(FREQ,X,Y,Z) and IMAG(FREQ,X,Y,Z)
%                     (This is the style in READ_DDF_CMPLX)
%
%   Note that this code can be used on gzip'd files.  These files are 
%   unzipped and reszipped during the reading process, using a shell
%   command.  This can be problematic if you are running on a Windows 
%   PC and the call SYSTEM('GUNZIP') doesn't mean anything to your PC.
%
%   Note that in general, this will NOT work on v3 files.
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
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/read_ddf_image.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%



if (nargin < 1)
  error('Not enough input arguments --- must at least give rootname');
end

if (nargin < 2) displayflag   = 1; end;
if (nargin < 3) complexformat = 0; end;

rootname = char(rootname);

unzipped_flag = 0;

if (displayflag == 1)
  disp(sprintf('\n  Reading DDF file:    %s.ddf',rootname));
end

output.ddf = readddf_file(rootname);

datafilename = strcat(rootname,'.cmplx');

datafile = fopen(datafilename,'r','b');
if (datafile < 0) 
  unzipped_flag = 1;
  if (exist(sprintf('%s.gz', datafilename)))
    unix(sprintf('gunzip %s.gz',datafilename));
    datafile = fopen(strcat(rootname,'.cmplx'),'r','b');  
  end
end
 
if (datafile < 0) 
  error('Could not find data file "%s" or "%s.gz"',datafilename, ...
        datafilename);
end

cmplxinfo = dir(datafilename);
if (cmplxinfo.bytes ~= prod(output.ddf.npix)*output.ddf.specpoints*2*4) 
  error('Expected "%s" to be %d bytes; found %d bytes on disk', ...
        datafilename, prod(output.ddf.npix)*output.ddf.specpoints*2*4, ...
        cmplxinfo.bytes);
end

if (displayflag == 1)
  disp(sprintf('  Reading data file:   %s',datafilename));
end

if (complexformat == 0)
  temp_data = fread(datafile, prod(output.ddf.npix)*output.ddf.specpoints*2, ...
        'float');
  if (numel(temp_data) ~= prod(output.ddf.npix)*output.ddf.specpoints*2)
    error('Only able to read %d floats, expected %d', ...
           numel(temp_data), ...
           prod(output.ddf.npix)*output.ddf.specpoints*2);    
  end
  output.img = zeros(length(temp_data)/2,1);
  output.img = complex(temp_data(1:2:length(temp_data)), ...
                       temp_data(2:2:length(temp_data)));
  output.img = reshape(output.img, ...
                       [output.ddf.specpoints output.ddf.npix]);
elseif (complexformat == 1)
  for (z = 1:output.ddf.npix(3))
    for (y = 1:output.ddf.npix(2))
      for (x = 1:output.ddf.npix(1))
        temp_data = fread(datafile, ...
                          output.ddf.specpoints*2, 'float');
        if (numel(temp_data) ~= 2*output.ddf.specpoints)
          error('cmplx file seems too short when reading (%d,%d,%d)', ... 
                 x, y, z);
        end

        if (max(abs(temp_data)) > eps)
          output.img(x,y,z).real = temp_data(1:2:end);
          output.img(x,y,z).imaginary = temp_data(2:2:end);    
        else
          output.img(x,y,z).real = [];
          output.img(x,y,z).imaginary = [];    
        end
      end
    end
  end
elseif (complexformat == 2)
  temp_data = fread(datafile, prod(output.ddf.npix)*output.ddf.specpoints*2, ...
                    'float');
  if (numel(temp_data) ~= prod(output.ddf.npix)*output.ddf.specpoints*2)
    error('Only able to read %d floats, expected %d', ...
           numel(temp_data), ...
           prod(output.ddf.npix)*output.ddf.specpoints*2);    
  end
  temp_complex = zeros(length(temp_data)/2,1);
  temp_complex = complex(temp_data(1:2:length(temp_data)), ...
                       temp_data(2:2:length(temp_data)));
  temp_complex = reshape(temp_complex, ...
                       [output.ddf.specpoints output.ddf.npix]);
  output.real = real(temp_complex);
  output.imag = imag(temp_complex);
end

if (unzipped_flag == 1)
  unix(sprintf('gzip %s', datafilename));
end

fclose(datafile);
