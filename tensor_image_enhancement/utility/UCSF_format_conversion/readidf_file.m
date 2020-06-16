function idf = readidf_file(rootname)
%
%   READIDF_FILE Reads an IDF file (but not the image file)
%
%   READIDF_FILE(ROOTNAME) reads in ROOTNAME.IDF and stores the
%   data as a Matlab structure.  This structure can be viewed easily
%   in Matlab, and can be passed to WRITEIDF_FILE or
%   WRITE_IDF_IMAGE for automated writing to disk.  Works for v3 or
%   v5 style IDFs.
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
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/readidf_file.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%



IDFname = strcat(rootname, '.idf');
fid = fopen(IDFname,'r','b');

if (fid < 0) 
  disp(sprintf('\n-- Error!  File %s not found! --\n',IDFname));
  return;
end

idf.filename = strcat(IDFname);

% 'IMAGE DESCRIPTOR FILE version'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('IMAGE DESCRIPTOR FILE version',2));
idf.version = sscanf(thisline(nextindex:end), '%d');

% 'studyid:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('studyid:',2));
idf.studyid = sscanf(thisline(nextindex:end), '%c');

% 'study #:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('study #:',2));
idf.studynum = sscanf(thisline(nextindex:end), '%c');

% 'series #:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('series #:',2));
idf.series = sscanf(thisline(nextindex:end), '%d');

% 'position:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('position: ',2));
idf.position = cleanupcomment(sscanf(thisline(nextindex:end), '%c'));

% 'coil:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('coil: ',2));
idf.coil = cleanupcomment(sscanf(thisline(nextindex:end), '%c'));

% 'orientation:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('orientation: ',2));
idf.orientation = sscanf(thisline(nextindex:end), '%d');

% 'echo/time/met index:'  ... and ... 'value:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('echo/time/met index:',2));
[idf.index, count, errmsg, nextindex] = ...
    sscanf(thisline(nextindex:end), '%d');
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', ...
	   size('echo/time/met index:           value:',2));
[idf.value, count, errmsg, nextindex] = ...
    sscanf(thisline(nextindex:end), '%f');

% 'rootname:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('rootname: ',2));
idf.rootname = cleanupcomment(sscanf(thisline(nextindex:end), ...
				     '%c'));

% 'comment:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('comment: ',2));
idf.comment = cleanupcomment(sscanf(thisline(nextindex:end), ...       
				    '%c'));

% 'filetype:'  ... and ...  'entry/pixel:' 

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('filetype:',2));
idf.filetype = sscanf(thisline(nextindex:end), '%d');
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('filetype:         entry/pixel:',2));
idf.entriesperpixel = ...
    sscanf(thisline(nextindex:end), '%d');

% dimensions

for (i = 1:3)

  thisline = fgetl(fid);

  thisline = fgetl(fid);

  [junk, count, errmsg, nextindex] = ...
      sscanf(thisline, '%c', ...
	     size('npix:',2));
  idf.npix(i) = sscanf(thisline(nextindex:end), '%d');
  
  [junk, count, errmsg, nextindex] = ...
      sscanf(thisline, '%c', ...
	     size('npix:         fov(mm):',2));
  idf.fov(i) = sscanf(thisline(nextindex:end), '%f');
  
  [junk, count, errmsg, nextindex] = ...
      sscanf(thisline, '%c', ...
	     size('npix:         fov(mm):          center(mm):',2));
  idf.center(i) = sscanf(thisline(nextindex:end), '%f');
  
  [junk, count, errmsg, nextindex] = ...
      sscanf(thisline, '%c', ...
	     size(['npix:   256   fov(mm):  240.00  center(mm):' ...
		   '    0.00  pixelsize(mm): '],2)); 
  idf.pixelsize(i) = sscanf(thisline(nextindex:end), '%f');
  
end

% 'slice thickness (mm):'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('slice thickness (mm):',2));
idf.slicethickness = sscanf(thisline(nextindex:end), '%f');

% 'minimum:'  ... and ...  'maximum:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('minimum:',2));
idf.minimum = sscanf(thisline(nextindex:end), '%f');
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('minimum:                  maximum:',2));
idf.maximum = ...
    sscanf(thisline(nextindex:end), '%f');

% 'scale:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('scale:',2));
idf.scale = sscanf(thisline(nextindex:end), '%f');

% 'first slice read:'

thisline = fgetl(fid);
[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size('first slice read:',2));
idf.firstread = sscanf(thisline(nextindex:end), '%d');

[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size(['first slice read:        last slice' ...
		    ' read:'],2));
idf.lastread = sscanf(thisline(nextindex:end), '%d');

[junk, count, errmsg, nextindex] = ...
    sscanf(thisline, '%c', size(['first slice read:        last slice' ...
		    ' read:        sliceskip:'],2));
idf.skip = sscanf(thisline(nextindex:end), '%d');

% The following is only read for v5 data

if (idf.version == 5) 

  thisline = fgetl(fid);

  if (isstr(thisline)) 
  
    thisline = fgetl(fid);
    [junk, count, errmsg, nextindex] = ...
        sscanf(thisline, '%c', size('center:',2));
    idf.LPScenter(1:3) = sscanf(thisline(nextindex:end), '%f %f %f');
  
    thisline = fgetl(fid);
    [junk, count, errmsg, nextindex] = ...
        sscanf(thisline, '%c', size('toplc:',2));
    idf.toplc(1:3) = sscanf(thisline(nextindex:end), '%f %f %f');
  
    thisline = fgetl(fid);
    [junk, count, errmsg, nextindex] = ...
        sscanf(thisline, '%c', size('dcos1:',2));
    idf.dcos(1,1:3) = sscanf(thisline(nextindex:end), '%f %f %f');
  
    thisline = fgetl(fid);
    [junk, count, errmsg, nextindex] = ...
        sscanf(thisline, '%c', size('dcos2:',2));
    idf.dcos(2,1:3) = sscanf(thisline(nextindex:end), '%f %f %f');

    thisline = fgetl(fid);
    [junk, count, errmsg, nextindex] = ...
        sscanf(thisline, '%c', size('dcos3:',2));
    idf.dcos(3,1:3) = sscanf(thisline(nextindex:end), '%f %f %f');

  end

end
  
fclose(fid);
