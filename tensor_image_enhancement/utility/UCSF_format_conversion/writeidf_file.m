function writeidf_file(rootname, idf)
%
%   WRITEIDF_FILE Writes an IDF file (but not the image file)
%
%   WRITEIDF_FILE(ROOTNAME,IDF) writes ROOTNAME.IDF based on the
%   information in the IDF Matlab structure.  Works for v5 data only.
%
%   It is recommended that IDF be read in by READIDF_FILE or
%   READ_IDF_IMAGE and modified by hand in Matlab, rather than created
%   completely manually.
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
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/writeidf_file.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%



IDFname = strcat(rootname, '.idf');
fid = fopen(IDFname,'w','b');

%idf.version = 5;

if (fid < 0) 
  disp(sprintf('\n-- Error!  File %s not found! --\n',IDFname));
  return;
end

% Strip everything before the last slash

last_slash = max(find(rootname == '/'));
if (last_slash == length(rootname))
  error('writeidf_file: specified "%s" cannot be a rootname', rootname);
end

if (~isempty(last_slash))
  idf.rootname = rootname(last_slash+1:end);
else
  idf.rootname = rootname;
end

%idf.filename = rootname;
%idf.rootname = idf.filename;
%for (i = length(rootname):-1:1)
%  if (rootname(i) == '/')      
%    idf.rootname = rootname(i+1:end);
%    break;
%  end
%end

% 'IMAGE DESCRIPTOR FILE version'

fprintf(fid, 'IMAGE DESCRIPTOR FILE version %d\n', idf.version);
if (~isempty(idf.studyid))
  fprintf(fid, 'studyid: %s\n', idf.studyid);
else
  fprintf(fid, 'studyid:\n');
end
if (~isempty(idf.studyid))
  fprintf(fid, 'study #:%s\n', idf.studynum);
else
  fprintf(fid, 'study #:\n');
end
fprintf(fid, 'series #:%8d\n', idf.series);
fprintf(fid, 'position: %s\n', idf.position);
fprintf(fid, 'coil: %s\n', idf.coil);
fprintf(fid, 'orientation:%4d', idf.orientation);
if (idf.orientation == 13)
  fprintf(fid, '     axial normal\n');
elseif (idf.orientation == 12)
  fprintf(fid, '     coronal normal\n');
elseif (idf.orientation == 11)
  fprintf(fid, '     sagittal normal\n');
elseif (idf.orientation == 9)
  fprintf(fid, '     axial oblique\n');   % oblique plane added - JML
else
  fprintf(fid, '     other\n');
end
fprintf(fid,'echo/time/met index:%6d     value:%11.2f\n', ...
	idf.index, idf.value);
fprintf(fid,'rootname: %s\n', idf.rootname);
fprintf(fid,'comment: %s\n', idf.comment);
fprintf(fid,'filetype:%4d', idf.filetype);
fprintf(fid,'     entry/pixel:%3d', idf.entriesperpixel);
fprintf(fid,'     DICOM format images\n');

fprintf(fid, 'dimension:%3d     columns     itype:%3d\n', 1, 1);
fprintf(fid, 'npix:%6d   fov(mm):%8.2f  center(mm):%8.2f  pixelsize(mm):%11.5f\n', ...
	idf.npix(1), idf.fov(1), idf.center(1), idf.pixelsize(1));
fprintf(fid, 'dimension:%3d     rows        itype:%3d\n', 2, 2);
fprintf(fid, 'npix:%6d   fov(mm):%8.2f  center(mm):%8.2f  pixelsize(mm):%11.5f\n', ...
	idf.npix(2), idf.fov(2), idf.center(2), idf.pixelsize(2));
fprintf(fid, 'dimension:%3d     slices      itype:%3d\n', 3, 3);
fprintf(fid, 'npix:%6d   fov(mm):%8.2f  center(mm):%8.2f  pixelsize(mm):%11.5f\n', ...
	idf.npix(3), idf.fov(3), idf.center(3), idf.pixelsize(3));


fprintf(fid,'slice thickness (mm):%15.5f\n', idf.slicethickness);

if (abs(idf.minimum) > 1e5) 
  fprintf(fid,'minimum:  %11.4g     maximum:  %11.4g\n', idf.minimum, idf.maximum);
else
  fprintf(fid,'minimum: %8.2f         maximum: %12.4g\n', idf.minimum, idf.maximum);
end

if (abs(idf.minimum) > 1e6) 
  fprintf(fid,'scale:%13.4g\n', idf.scale);
else
  fprintf(fid,'scale:%10.3f\n', idf.scale);
end

fprintf(fid, 'first slice read:%5d   last slice read:%5d   sliceskip:%5d\n', ...
	idf.firstread, idf.lastread, idf.skip);

% The old version forced it to be axial SI oriented.  No more!

%idf.LPScenter = [idf.center(1), idf.center(2), -idf.center(3)];
%idf.toplc = [idf.center(1) - (idf.fov(1)/2) + idf.pixelsize(1)/2 ...
%	     idf.center(2) - (idf.fov(2)/2) + idf.pixelsize(2)/2 ...
%	     -idf.center(3) + (idf.fov(3)/2) - idf.pixelsize(3)/2];
%idf.dcos1 = [1 0 0];
%idf.dcos2 = [0 1 0];
%idf.dcos3 = [0 0 -1];


if (idf.version == 5) 
  
  fprintf(fid, 'LOCATION DATA IN LPS COORDINATES\n');  
  fprintf(fid, 'center: %14.5f%14.5f%14.5f\n', idf.LPScenter);
  fprintf(fid, 'toplc:  %14.5f%14.5f%14.5f\n', idf.toplc);
  fprintf(fid, 'dcos1:  %14.5f%14.5f%14.5f\n', idf.dcos(1,:));  
  fprintf(fid, 'dcos2:  %14.5f%14.5f%14.5f\n', idf.dcos(2,:));  
  fprintf(fid, 'dcos3:  %14.5f%14.5f%14.5f\n', idf.dcos(3,:));    

end
  
fclose(fid);
