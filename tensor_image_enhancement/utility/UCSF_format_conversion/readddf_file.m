function ddf = readddf_file(rootname)
%
%   READDDF_FILE Reads a DDF file (but not the complex file)
%
%   READDDF_FILE(ROOTNAME) reads in ROOTNAME.DDF and stores the
%   data as a Matlab structure.  This structure can be viewed easily
%   in Matlab, and can be passed to WRITEDDF_FILE or 
%   WRITE_DDF_IMAGE for automated writing to disk.  Works for v3, v5 or 
%   v6 style DDFs.
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
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/readddf_file.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%


if (nargin == 0)
  error('You must specify the rootname as the input argument');
end

DDFname = strcat(rootname, '.ddf');
fid = fopen(DDFname,'r','b');

if (fid < 0) 
  error('Unable to open specified file: %s', DDFname);
end

ddf.filename = DDFname;

v5a_field_defs = ...
    {'version',                       '%f',  inf, 'version'         , []  
     'study id:',                     '%c',  inf, 'studyid'         , []
     'study # :',                     '%c',  inf, 'studynum'        , []
     'position:',                     '%c',  inf, 'position'        , []
     'coil:',                         '%c',  inf, 'coil'            , []
     'series #:'                      '%d',  inf, 'series'          , []
     'orientation:'                   '%d',  inf, 'orientation'     , []
     'echo/time/met index:',          '%d',  inf, 'index'           , []
     'value:',                        '%f',  inf, 'value'           , []
     'name:',                         '%c',  inf, 'rootname'        , []
     'comment:',                      '%c',  inf, 'comment'         , []
     'sequence name:',                '%c',  inf, 'sequence'        , []
     'localization type:',            '%d',  inf, 'localizationtype', []
     ':',                             '%f',  inf, 'specfrequency'   , []
     'sweepwidth(Hz):',               '%f',  inf, 'sweepwidth'      , []
     'dwelltime:',                    '%f',  inf, 'dwelltime'       , []
     'satpulse pos(Hz):',             '%f',  inf, 'satpulsepos'     , []
     'bandwidth(Hz):',                '%f',  inf, 'bandwidth'       , []
     'beginning of acquisition (ms):','%f',  inf, 'acqbeginning'    , []
     'gradient on time (ms):',        '%f',  inf, 'gradontime'      , []
     'nex:',                          '%f',  inf, 'nex'             , []
     'chop:',                         '%f',  inf, 'chop'            , []
     'even_sym:',                     '%f',  inf, 'even_sym'        , []
     'pe_order:',                     '%f',  inf, 'pe_order'        , []
     'datatype:',                     '%d',  inf, 'datatype'        , []
     'number of dimensions:',         '%d',  inf, 'numdim'          , []
     'npts:',                         '%d',  inf, 'specpoints'      , []
     'acqpts:',                       '%d',  inf, 'acqspec'         , []
     'npts:',                         '%d',  inf, 'npix'            , 1
     'fov:',                          '%f',  inf, 'fov'             , 1
     'acqpts:',                       '%d',  inf, 'acqpts'          , 1     
     'npts:',                         '%d',  inf, 'npix'            , 2
     'fov:',                          '%f',  inf, 'fov'             , 2
     'acqpts:',                       '%d',  inf, 'acqpts'          , 2
     'npts:',                         '%d',  inf, 'npix'            , 3
     'fov:',                          '%f',  inf, 'fov'             , 3     
     'acqpts:',                       '%d',  inf, 'acqpts'          , 3     
     'center:',                       '%f',[1,3], 'recon_center'    , []
     'toplc:',                        '%f',[1,3], 'recon_toplc'     , []
     'dcos:',                         '%f',[1,3], 'recon_dcos1'     , []
     'dcos:',                         '%f',[1,3], 'recon_dcos2'     , []
     'dcos:',                         '%f',[1,3], 'recon_dcos3'     , []
     'pcenter:',                      '%f',[1,3], 'encode_center'   , []
     'pfov:',                         '%f',[1,3], 'encode_fov'      , []
     'pmatrix:',                      '%d',[1,3], 'encode_matrix'   , []
     'reverse:',                      '%d',[1,3], 'encode_reverse'  , []
     'pdcos:',                        '%f',[1,3], 'encode_dcos1'    , []
     'pdcos:',                        '%f',[1,3], 'encode_dcos2'    , []
     'pdcos:',                        '%f',[1,3], 'encode_dcos3'    , []
     'bcenter:',                      '%f',[1,3], 'box_center'      , []
     'bsize:',                        '%f',[1,3], 'box_size'        , []
     'bdcos:',                        '%f',[1,3], 'box_dcos1'       , []
     'bdcos:',                        '%f',[1,3], 'box_dcos2'       , []
     'bdcos:',                        '%f',[1,3], 'box_dcos3'       , []
     'refacq:',                       '%c',  inf, 'refacq'          , []
     'refrecon:',                     '%c',  inf, 'refrecon'        , []
    };

v5b_field_defs = ...
    {'version',                       '%f',  inf, 'version'         , []  
     'study id:',                     '%c',  inf, 'studyid'         , []
     'study # :',                     '%c',  inf, 'studynum'        , []
     'position:',                     '%c',  inf, 'position'        , []
     'coil:',                         '%c',  inf, 'coil'            , []
     'series #:'                      '%d',  inf, 'series'          , []
     'orientation:'                   '%d',  inf, 'orientation'     , []
     'echo/time/met index:',          '%d',  inf, 'index'           , []
     'value:',                        '%f',  inf, 'value'           , []
     'name:',                         '%c',  inf, 'rootname'        , []
     'comment:',                      '%c',  inf, 'comment'         , []
     'sequence name:',                '%c',  inf, 'sequence'        , []
     'localization type:',            '%d',  inf, 'localizationtype', []
     ':'             ,                '%f',  inf, 'specfrequency'   , []
     'sweepwidth(Hz):',               '%f',  inf, 'sweepwidth'      , []
     'dwelltime:',                    '%f',  inf, 'dwelltime'       , []
     'centfreq pos(Hz):',             '%f',  inf, 'centfreq'        , []
     'pulse type:',                   '%d',  inf, 'pulsetype'        , []
     'beginning of acquisition (ms):','%f',  inf, 'acqbeginning'    , []
     'gradient on time (ms):',        '%f',  inf, 'gradontime'      , []
     'nex:',                          '%f',  inf, 'nex'             , []
     'chop:',                         '%f',  inf, 'chop'            , []
     'even_sym:',                     '%f',  inf, 'even_sym'        , []
     'pe_order:',                     '%f',  inf, 'pe_order'        , []
     'datatype:',                     '%d',  inf, 'datatype'        , []
     'number of dimensions:',         '%d',  inf, 'numdim'          , []
     'npts:',                         '%d',  inf, 'specpoints'      , []
     'acqpts:',                       '%d',  inf, 'acqspec'         , []
     'npts:',                         '%d',  inf, 'npix'            , 1
     'fov:',                          '%f',  inf, 'fov'             , 1
     'acqpts:',                       '%d',  inf, 'acqpts'          , 1     
     'npts:',                         '%d',  inf, 'npix'            , 2
     'fov:',                          '%f',  inf, 'fov'             , 2
     'acqpts:',                       '%d',  inf, 'acqpts'          , 2
     'npts:',                         '%d',  inf, 'npix'            , 3
     'fov:',                          '%f',  inf, 'fov'             , 3     
     'acqpts:',                       '%d',  inf, 'acqpts'          , 3     
     'center:',                       '%f',[1,3], 'recon_center'    , []
     'toplc:',                        '%f',[1,3], 'recon_toplc'     , []
     'dcos:',                         '%f',[1,3], 'recon_dcos1'     , []
     'dcos:',                         '%f',[1,3], 'recon_dcos2'     , []
     'dcos:',                         '%f',[1,3], 'recon_dcos3'     , []
     'pcenter:',                      '%f',[1,3], 'encode_center'   , []
     'pfov:',                         '%f',[1,3], 'encode_fov'      , []
     'pmatrix:',                      '%d',[1,3], 'encode_matrix'   , []
     'reverse:',                      '%d',[1,3], 'encode_reverse'  , []
     'pdcos:',                        '%f',[1,3], 'encode_dcos1'    , []
     'pdcos:',                        '%f',[1,3], 'encode_dcos2'    , []
     'pdcos:',                        '%f',[1,3], 'encode_dcos3'    , []
     'bcenter:',                      '%f',[1,3], 'box_center'      , []
     'bsize:',                        '%f',[1,3], 'box_size'        , []
     'bdcos:',                        '%f',[1,3], 'box_dcos1'       , []
     'bdcos:',                        '%f',[1,3], 'box_dcos2'       , []
     'bdcos:',                        '%f',[1,3], 'box_dcos3'       , []
     'refacq:',                       '%c',  inf, 'refacq'          , []
     'refrecon:',                     '%c',  inf, 'refrecon'        , []
    };

v6_field_defs = ...
    {'version:',                      '%f',  inf, 'version'         , []
     'object type:',                  '%c',  inf, 'object_type'     , []
     'patient id:',                   '%c',  inf, 'patient_id'      , []
     'patient name:',                 '%c',  inf, 'patient_name'    , []
     'patient code:',                 '%c',  inf, 'patient_code'    , []
     'date of birth:',                '%c',  inf, 'dob'             , []
     'sex:',                          '%c',  inf, 'sex'             , []
     'study id:',                     '%c',  inf, 'studyid'         , []
     'study code:',                   '%c',  inf, 'study_code'      , []
     'study date:',                   '%c',  inf, 'study_date'      , []
     'accession number:',             '%c',  inf, 'accession_number', []
     'name:',                         '%c',  inf, 'root_name'       , []
     'series number:',                '%d',  inf, 'series'          , []
     'series description:'            '%c',  inf, 'series_description', []
     'comment:'                       '%c',  inf, 'comment'         , []
     'patient entry:',                '%c',  inf, 'patient_entry'   , []
     'patient position:',             '%c',  inf, 'patient_position', []
     'orientation:',                  '%c',  inf, 'orientation'     , []
     'data type:',                    '%c',  inf, 'data_type'       , []
     'number of components:',         '%d',  inf, 'num_components'  , []
     'source description:',           '%c',  inf, 'source_description', []
     'number of dimensions:',         '%d',  inf, 'numdim'          , []
     'type:',                         '%c',  inf, 'dimension_type'  , 1
     'npoints:',                      '%d',  inf, 'specpoints'      , []
     'type:',                         '%c',  inf, 'dimension_type'  , 2
     'npoints:',                      '%d',  inf, 'npix'            , 1
     'pixel spacing(mm):',            '%f',  inf, 'pixel_spacing'   , 1     
     'type:',                         '%c',  inf, 'dimension_type'  , 3
     'npoints:',                      '%d',  inf, 'npix'            , 2
     'pixel spacing(mm):',            '%f',  inf, 'pixel_spacing'   , 2     
     'type:',                         '%c',  inf, 'dimension_type'  , 4
     'npoints:',                      '%d',  inf, 'npix'            , 3
     'pixel spacing(mm):',            '%f',  inf, 'pixel_spacing'   , 3 
     'center(lps, mm):',              '%f',[1,3], 'center'          , []
     'toplc(lps, mm):',               '%f',[1,3], 'toplc'           , []
     'dcos0:',                        '%f',[1,3], 'dcos0'           , []
     'dcos1:',                        '%f',[1,3], 'dcos1'           , []
     'dcos2:',                        '%f',[1,3], 'dcos2'           , []
     'coil name:',                    '%c',  inf, 'coil_name'       , []
     'slice gap(mm):',                '%f',  inf, 'slice_gap'       , []
     'echo time(ms):',                '%f',  inf, 'TE'              , []
     'repetition time(ms):',          '%f',  inf, 'TR'              , []     
     'inversion time(ms):',           '%f',  inf, 'TI'              , []
     'flip angle:',                   '%f',  inf, 'flip_angle'      , []
     'pulse sequence name:',          '%c',  inf, 'sequence'        , []
     'transmitter frequency(MHz):',   '%f',  inf, 'transmit_freq'   , []
     'isotope:',                      '%c',  inf, 'isotope'         , []
     'field strength(T):',            '%f',  inf, 'field_strength'  , []
     'number of sat bands:',          '%d',  inf, 'num_sat_bands'   , []
     'localization type:',            '%c',  inf, 'loc_type'        , []
     'center frequency(MHz):',        '%f',  inf, 'centfreq'        , []
     'ppm reference:',                '%f',  inf, 'ppm_ref'         , []
     'sweepwidth(Hz):',               '%f',  inf, 'sweepwidth'      , []
     'dwelltime(ms):',                '%f',  inf, 'dwelltime'       , []
     'frequency offset(Hz):',         '%f',  inf, 'freq_offset'     , []
     'centered on water:',            '%c',  inf, 'center_on_water' , []
     'suppression technique:',        '%c',  inf, 'suppresion_tech' , []
     'residual water:',               '%c',  inf, 'residual_water'  , []
     'number of acquisitions:',       '%d',  inf, 'num_acquisitions', []
     'chop:',                         '%c',  inf, 'chop'            , []
     'even symmetry:',                '%c',  inf, 'even_sym'        , []
     'data reordered:',               '%c',  inf, 'data_reordered'  , []
     'acq. toplc(lps, mm):',          '%f',[1,3], 'acq_toplc'       , []
     'acq. spacing(mm):',             '%f',[1,3], 'acq_spacing'     , []
     'acq. number of data points:',   '%d',  inf, 'acq_n_data_points'   , []
     'acq. number of points:',        '%d',[1,3], 'acq_n_points'    , []
     'acq. dcos1:',                   '%f',[1,3], 'acq_dcos1'       , []
     'acq. dcos2:',                   '%f',[1,3], 'acq_dcos2'       , []
     'acq. dcos3:',                   '%f',[1,3], 'acq_dcos3'       , []
     'selection center(lps, mm):',    '%f',[1,3], 'box_center'      , []     
     'selection size(mm):',           '%f',[1,3], 'box_size'        , []
     'selection dcos1:',              '%f',[1,3], 'box_dcos1'       , []
     'selection dcos2:',              '%f',[1,3], 'box_dcos2'       , []
     'selection dcos3:',              '%f',[1,3], 'box_dcos3'       , []
    };

v6p1_field_defs = cat(1, v6_field_defs, ...
    {
     'reordered toplc(lps, mm):',     '%f',[1,3], 'reordered_toplc'      , []
     'reordered center(lps, mm):',    '%f',[1,3], 'reordered_center'     , []
     'reordered spacing(mm):',        '%f',[1,3], 'reordered_spacing'    , []
     'reordered number of points:',   '%d',[1,3], 'reordered_n_points'   , []
     'reordered dcos1:',              '%f',[1,3], 'reordered_dcos1'      , []
     'reordered dcos2:',              '%f',[1,3], 'reordered_dcos2'      , []
     'reordered dcos3:',              '%f',[1,3], 'reordered_dcos3'      , []
    });


entire_file = fscanf(fid,'%c');
fclose(fid);

% First thing to do is find the version number

version = sscanf(entire_file(strfind(entire_file, 'version')+8:end), '%f');

if (version == 3)  % Convert to v5   
  temp_input_name = sprintf('/tmp/read_ddf_input_%010d.txt',sum(clock));
  temp_ddf_name   = sprintf('/tmp/read_ddf_v5_%010d',sum(clock));
  temp_input = fopen(temp_input_name,'w');
  fprintf(temp_input, '%s.ddf\n%s\n', rootname, temp_ddf_name);
  fclose(temp_input);
  system(sprintf('convert_ddf_v5 < %s > /dev/null', temp_input_name));
  delete(temp_input_name);  
  ddf = readddf_file(temp_ddf_name);
  ddf.filename = strcat(DDFname);
  delete(strcat(temp_ddf_name,'.ddf')); 
  return;
elseif (abs(version - 5) < 0.5) 
  if (~isempty(strfind(entire_file,'satpulse')))
    field_defs = v5a_field_defs;
  else
    field_defs = v5b_field_defs;
  end
elseif (version == 6)
  field_defs = v6_field_defs;
elseif (version == 6.1)
  field_defs = v6p1_field_defs;
else
  error('Sorry, ddf version %g is not supported\n', version);
end

number_of_fields = size(field_defs,1);

location_end = 1;

% =======================================================================
% Begin reading the actual data
% =======================================================================

for (fieldnum = 1:number_of_fields)
    
  % The data is located between this field name and the next field name:
  % the fields in 'field_defs' are given in order, so only have to look
  % in a small region to find the field we want.  This is important
  % since the DDF has a lot of other information about the
  % processing at the end of the DDF ... consequently, there may be
  % multiple instances of fields ... for example it may specify
  % 'dcos' 3 or 4 times within a given DDF file, but we only want
  % to read in a specific instance of this.
  
  location_start = ...
      min(strfind(entire_file(location_end:end), ...
              field_defs{fieldnum,1})) + length(field_defs{fieldnum,1}) ...
      + location_end;
  
  if (isempty(location_start))
    error('Unable to find field "%s" in file "%s"', ...
          field_defs{fieldnum,1}(1:end-1), DDFname);
  end
  
  if (fieldnum ~= number_of_fields)  
    location_end = min(strfind(entire_file(location_start:end), ...
                               field_defs{fieldnum+1,1}))-1+location_start;
  else
    location_end = length(entire_file)+1;
  end
  
  % Read in the data
  
  field_name  = field_defs{fieldnum,4};  
  field_value = sscanf(entire_file(location_start:location_end-1), ...
                       field_defs{fieldnum,2},                   ...
                       field_defs{fieldnum,3});
  
  % Clean up the data, in the case of strings
  
  if (isstr(field_value))
    field_value = trimWhitespace(field_value);  
    if (isempty(field_value))
      field_value = '';
    end
  end
  
  % Assign the data into the structure
  
  array_index = field_defs{fieldnum,5};
  
  if (isempty(array_index))
    ddf = setfield(ddf, field_name, field_value);
  else
    if (isfield(ddf, field_name))
      temp_value = getfield(ddf, field_name);
    else
      temp_value = [];    
    end
    if (isstr(field_value)) 
      temp_value{array_index} = field_value;
    else
      temp_value(array_index) = field_value;      
    end
    ddf = setfield(ddf, field_name, temp_value);
  end
  
end

if (version >= 6)

  % For v6 data, we read in the sat band information.  We weren't
  % able to do this beforehand, because the number of sat bands was
  % not determined until we've read in the bulk of the DDF file.
  
  for (satnum = 1:ddf.num_sat_bands)
    
    % Find the location of the line "sat band n thickness(mm):"
    
    location_start = ...
        min(strfind(entire_file, ...
                    sprintf('sat band %d thickness(mm):', satnum))) + ...
        length(sprintf('sat band %d thickness(mm):', satnum));
    location_end = ...
        min(strfind(entire_file, ...
                    sprintf('sat band %d orientation:',satnum)))-1;
    
    ddf.satband_thickness(satnum,1) = ...
        sscanf(entire_file(location_start:location_end), '%f',satnum);
    
    % Find the location of the line "sat band 1 thickness(mm):"
    
    location_start = ...
        min(strfind(entire_file, ...
                    sprintf('sat band %d orientation:',satnum))) + ...
        length(sprintf('sat band %d orientation:',satnum));
    location_end = ...
        min(strfind(entire_file, ...
                    sprintf('sat band %d position(lps, mm):',satnum)))-1;
    
    ddf.satband_orientation(satnum,:) = ...
        sscanf(entire_file(location_start:location_end), '%f', [1,3]);
    
    % Find the location of the line "sat band 1 position(lps, mm):"
    
    location_start = ...
        min(strfind(entire_file, ...
                    sprintf('sat band %d position(lps, mm):',satnum)))+...
        length(sprintf('sat band %d position(lps, mm):',satnum));
    location_end = length(entire_file);

    ddf.satband_position(satnum,:) = ...
        sscanf(entire_file(location_start:location_end), '%f', [1,3]);  
  
  end

  % For v6, we also read in all the other "stuff" that is tacked on
  % to the end of the DDF.  This typically has information about
  % all the processing that happened to be used to create this
  % file.  Thus, it has no predefined form, so we just grab
  % everything and dump it into a giant string so we can write it
  % back out or access it later, if required.
  
  separator = '===================================================';
  locations = strfind(entire_file, separator)+length(separator)+1;
  
  if (length(locations) >= 3)
    ddf.misc_info = entire_file(locations(3):end);
  else
    ddf.misc_info = [];
  end
end

% ddf.read_in_time = datestr(clock,0);

%ddf.npix = ddf.npts; % For backwards compatibility!

if (ddf.version >= 6)
  ddf.encode_fov = ddf.npix .* ddf.pixel_spacing;   % For backwards compatibility!
end

% ====================================

function input_string = trimWhitespace(input_string)

if (~isempty(input_string))
  first_newline = min(find(input_string == sprintf('\n')));
  if (~isempty(first_newline))
    input_string = input_string(1:first_newline-1);
  end
end

if (~isempty(input_string))
  while (input_string(end) == sprintf('\n') | input_string(end) == ' ')
    input_string = input_string(1:end-1);
    if (isempty(input_string)) 
      break;
    end
  end
end

if (~isempty(input_string))
  while (input_string(1) == sprintf('\n') | input_string(1) == ' ')
    input_string = input_string(2:end);
    if (isempty(input_string)) 
      break;
    end
  end
end

