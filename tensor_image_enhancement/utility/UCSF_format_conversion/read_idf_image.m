function output = read_idf_image(rootname, displayflag, byte_order)
%
%   READ_IDF_IMAGE Reads an IDF file and associated image file
%
%   IMAGE_STRUCTURE = READ_IDF_IMAGE(ROOTNAME) reads two files for this
%   image.  First, it reads 'ROOTNAME.idf' which gives it the associated
%   header information.  This IDF can be either an MRSC v3 or v5 file.  
%   An image file is then read; these are .byt (unsigned 8-bit chars), 
%   .int2 (16-bit integers) or .real (32-bit floats).  The extension
%   that is attached to ROOTNAME is decided based on the 'FILETYPE' field
%   that is read out of the IDF file.  Note that the image file name
%   is based on the entered ROOTNAME plus the extension, not the rootname
%   that is given in the IDF file.  This is a departure from the way that
%   the IDL code MRSC_IMAGE works.  The return value is a structure with
%   two fields.  The first is IMAGE_STRUCTURE.IDF where IDF is itself a 
%   strcuture with fields corresponding to data entries found in the IDF
%   file (things like number of voxels and dimensions.)  The second
%   field in the IMAGE_STRUCTURE structure is the image itself, which is
%   just a 3D array of the actual voxel values.
%
%   IMAGE_STRUCTURE = READ_IDF_IMAGE(ROOTNAME, DISPLAYFLAG) does the same
%   thing, but if DISPLAYFLAG is a 0, then it doesn't print out information
%   about what file it's opening.  This is useful in those cases when you
%   are trying to have a nice clean uncluttered screen.
%
%   IMAGE_STRUCTURE = READ_IDF_IMAGE(ROOTNAME, DISPLAYFLAG, BYTE_ORDER) 
%   where BYTE_ORDER is 'n', 'b', or 'l', for native, big-endian,
%   or litte-endian byte-ordering, respectively (provide the letter
%   with the single quotes).  The default is 'b'.
%
%   Note that this code can be used on gzip'd files.  These files are 
%   unzipped and reszipped during the reading process, using a shell
%   command.  This can be problematic if you are running on a Windows 
%   PC and the call SYSTEM('GUNZIP') doesn't mean anything to your PC.
%
%   See also READIDF_FILE, WRITE_IDF_IMAGE, WRITEIDF_FILE.
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
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/read_idf_image.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%



if (nargin < 1)
  %help read_idf_image;
  error('Not enough input arguments --- must give at least rootname');
end

rootname = char(rootname);

if (nargin < 2)
  displayflag = 1;
end

if (isempty(displayflag))
  displayflag = 1;
end

if (nargin < 3)
  byte_order = 'b';
end

unzipped_flag = 0;

if (displayflag == 1)
  disp(sprintf('\n  Reading IDF file:    %s.idf',rootname));
end

output.idf = readidf_file(rootname);

switch output.idf.filetype
 case 2
  extension = '.byt';
  datatype  = 'uchar';
  bytespervox = 1;
 case 3
  extension = '.int2';
  datatype  = 'int16';
  bytespervox = 2;
 case 7
  extension = '.real';
  datatype  = 'float';
  bytespervox = 4;
 case 17
  extension = '.cmplx';
  datatype  = 'float';
  bytespervox = 8;
 otherwise
  error('Unrecognized file type "%d"', ...
        output.idf.filetype);
end

% Note that we are getting rootname from whatever the user entered,
% and NOT reading it from the IDF file.  If you want it to read from
% the IDF file, change ROOTNAME in this next line to OUTPUT.IDF.ROOTNAME.

imagefilename = strcat(rootname,extension);

imagefile = fopen(imagefilename,'r', byte_order);

if (imagefile < 0) 
  unzipped_flag = 1;
  if (exist(sprintf('%s.gz', imagefilename)))
    unix(sprintf('gunzip %s.gz',imagefilename));
    imagefile = fopen(strcat(rootname,extension),'r','b');  
  end
end

if (imagefile < 0)
  error('Could not find image file "%s" or "%s.gz"', ...
        imagefilename, imagefilename);
end

if (displayflag == 1) 
  disp(sprintf('  Reading image file:  %s',imagefilename));
end

fileinfo = dir(imagefilename);
if (fileinfo.bytes ~= prod(output.idf.npix)*bytespervox)
  error('Expected "%s" to be %d bytes; found %d bytes on disk', ...
        imagefilename, prod(output.idf.npix)*bytespervox, fileinfo.bytes);
end

if (output.idf.filetype ~= 17)   
  output.img = reshape(fread(imagefile, ...
                             prod(output.idf.npix), ...
                             datatype), output.idf.npix);
else
  temp = fread(imagefile, 2*prod(output.idf.npix), ...
               datatype);
  output.img = complex(reshape(temp(1:2:end), output.idf.npix), ...
                       reshape(temp(2:2:end), output.idf.npix));
end

if (unzipped_flag == 1)
  unix(sprintf('gzip %s', imagefilename));
end

fclose(imagefile);
