function ddf = writeddf_file(outputname, ddf)
%   WRITEDDF_FILE Writes a DDF file (but not the complex file)
%
%   WRITEDDF_FILE(ROOTNAME,DDF) writes ROOTNAME.DDF based on the 
%   information in the DDF Matlab structure.  Works for v5 or v6 data.
%
%   It is recommended that DDF be read in by READDDF_FILE or 
%   READ_DDF_IMAGE and modified by hand in Matlab, rather than created 
%   completely manually.  
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
%   $URL: https://intrarad.ucsf.edu/svn/rad_software/surbeck/brain/libs/file_io/trunk/writeddf_file.m $
%   $Rev: 14674 $
%   $Author: jasonc@RADIOLOGY.UCSF.EDU $
%   $Date: 2009-08-25 16:29:52 -0700 (Tue, 25 Aug 2009) $
%



v5a_field_defs = ...
    {'version',                       '%s',  inf, 'version'         , []  
     'study id:',                     '%s',  inf, 'studyid'         , []
     'study # :',                     '%s',  inf, 'studynum'        , []
     'position:',                     '%s',  inf, 'position'        , []
     'coil:',                         '%s',  inf, 'coil'            , []
     'series #:'                      '%d',  inf, 'series'          , []
     'orientation:'                   '%d',  inf, 'orientation'     , []
     'echo/time/met index:',          '%d',  inf, 'index'           , []
     'value:',                        '%f',  inf, 'value'           , []
     'root name:',                    '%c',  inf, 'rootname'        , []
     'comment:',                      '%c',  inf, 'comment'         , []
     'sequence name:',                '%c',  inf, 'sequence'        , []
     'localization type:',            '%d',  inf, 'localizationtype', []
     'spec frequency(MHz):',          '%f',  inf, 'specfrequency'   , []
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
    {'version',                       '%s',  inf, 'version'         , []  
     'study id:',                     '%c',  inf, 'studyid'         , []
     'study # :',                     '%c',  inf, 'studynum'        , []
     'position:',                     '%c',  inf, 'position'        , []
     'coil:',                         '%c',  inf, 'coil'            , []
     'series #:'                      '%d',  inf, 'series'          , []
     'orientation:'                   '%d',  inf, 'orientation'     , []
     'echo/time/met index:',          '%d',  inf, 'index'           , []
     'value:',                        '%f',  inf, 'value'           , []
     'root name:',                    '%c',  inf, 'rootname'        , []
     'comment:',                      '%c',  inf, 'comment'         , []
     'sequence name:',                '%c',  inf, 'sequence'        , []
     'localization type:',            '%d',  inf, 'localizationtype', []
     'spec frequency(kHz):',          '%f',  inf, 'specfrequency'   , []
     'sweepwidth(Hz):',               '%f',  inf, 'sweepwidth'      , []
     'dwelltime:',                    '%f',  inf, 'dwelltime'       , []
     'centfreq pos(Hz):',             '%f',  inf, 'centfreq'        , []
     'pulse type:',                   '%d',  inf, 'pulsetype'       , []
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
    {'version:',                      '%s',  inf, 'version'         , []
     'object type:',                  '%s',  inf, 'object_type'     , []
     'patient id:',                   '%s',  inf, 'patient_id'      , []
     'patient name:',                 '%s',  inf, 'patient_name'    , []
     'patient code:',                 '%s',  inf, 'patient_code'    , []
     'date of birth:',                '%s',  inf, 'dob'             , []
     'sex:',                          '%s',  inf, 'sex'             , []
     'study id:',                     '%s',  inf, 'studyid'         , []
     'study code:',                   '%s',  inf, 'study_code'      , []
     'study date:',                   '%s',  inf, 'study_date'      , []
     'accession number:',             '%s',  inf, 'accession_number', []
     'root name:',                    '%s',  inf, 'root_name'       , []
     'series number:',                '%d',  inf, 'series'          , []
     'series description:'            '%s',  inf, 'series_description', []
     'comment:'                       '%s',  inf, 'comment'         , []
     'patient entry:',                '%s',  inf, 'patient_entry'   , []
     'patient position:',             '%s',  inf, 'patient_position', []
     'orientation:',                  '%s',  inf, 'orientation'     , []
     'data type:',                    '%s',  inf, 'data_type'       , []
     'number of components:',         '%d',  inf, 'num_components'  , []
     'source description:',           '%s',  inf, 'source_description', []
     'number of dimensions:',         '%d',  inf, 'numdim'          , []
     'type:',                         '%s',  inf, 'dimension_type'  , 1
     'npoints:',                      '%d',  inf, 'specpoints'      , []
     'type:',                         '%s',  inf, 'dimension_type'  , 2
     'npoints:',                      '%d',  inf, 'npix'            , 1
     'pixel spacing(mm):',            '%f',  inf, 'pixel_spacing'   , 1     
     'type:',                         '%s',  inf, 'dimension_type'  , 3
     'npoints:',                      '%d',  inf, 'npix'            , 2
     'pixel spacing(mm):',            '%f',  inf, 'pixel_spacing'   , 2     
     'type:',                         '%s',  inf, 'dimension_type'  , 4
     'npoints:',                      '%d',  inf, 'npix'            , 3
     'pixel spacing(mm):',            '%f',  inf, 'pixel_spacing'   , 3 
     'center(lps, mm):',     '%14.5f%14.5f%14.5f',[1,3], 'center'   , []
     'toplc(lps, mm): ',     '%14.5f%14.5f%14.5f',[1,3], 'toplc'    , []
     'dcos0:',               '%14.5f%14.5f%14.5f',[1,3], 'dcos0'    , []
     'dcos1:',               '%14.5f%14.5f%14.5f',[1,3], 'dcos1'    , []
     'dcos2:',               '%14.5f%14.5f%14.5f',[1,3], 'dcos2'    , []
     'coil name:',                    '%s',  inf, 'coil_name'       , []
     'slice gap(mm):',                '%.6f',  inf, 'slice_gap'       , []
     'echo time(ms):',                '%.6f',  inf, 'TE'              , []
     'repetition time(ms):',          '%.6f',  inf, 'TR'              , []     
     'inversion time(ms):',           '%.6f',  inf, 'TI'              , []
     'flip angle:',                   '%.6f',  inf, 'flip_angle'      , []
     'pulse sequence name:',          '%s',  inf, 'sequence'        , []
     'transmitter frequency(MHz):',   '%.6f',  inf, 'transmit_freq'   , []
     'isotope:',                      '%s',  inf, 'isotope'         , []
     'field strength(T):',            '%.6f',  inf, 'field_strength'  , []
     'number of sat bands:',          '%d',  inf, 'num_sat_bands'   , []
     'localization type:',            '%s',  inf, 'loc_type'        , []
     'center frequency(MHz):',        '%.6f',  inf, 'centfreq'        , []
     'ppm reference:',                '%.6f',  inf, 'ppm_ref'         , []
     'sweepwidth(Hz):',               '%.6f',  inf, 'sweepwidth'      , []
     'dwelltime(ms):',                '%.6f',  inf, 'dwelltime'       , []
     'frequency offset(Hz):',         '%.6f',  inf, 'freq_offset'     , []
     'centered on water:',            '%s',  inf, 'center_on_water' , []
     'suppression technique:',        '%s',  inf, 'suppresion_tech' , []
     'residual water:',               '%s',  inf, 'residual_water'  , []
     'number of acquisitions:',       '%d',  inf, 'num_acquisitions', []
     'chop:',                         '%s',  inf, 'chop'            , []
     'even symmetry:',                '%s',  inf, 'even_sym'        , []
     'data reordered:',               '%s',  inf, 'data_reordered'  , []
     'acq. toplc(lps, mm):',  '%14.5f%14.5f%14.5f',[1,3], 'acq_toplc', []
     'acq. spacing(mm):   ',  '%14.5f%14.5f%14.5f',[1,3], 'acq_spacing', []
     'acq. number of data points:', '%d', inf, 'acq_n_data_points'   , []
     'acq. number of points:',  '%d %d %d',[1,3], 'acq_n_points'    , []
     'acq. dcos1:',      '%14.5f%14.5f%14.5f',[1,3], 'acq_dcos1'       , []
     'acq. dcos2:',      '%14.5f%14.5f%14.5f',[1,3], 'acq_dcos2'       , []
     'acq. dcos3:',      '%14.5f%14.5f%14.5f',[1,3], 'acq_dcos3'       , []
     'selection center(lps, mm):','%14.5f%14.5f%14.5f',[1,3],'box_center',[]
     'selection size(mm):       ','%14.5f%14.5f%14.5f',[1,3], 'box_size'    ,[]
     'selection dcos1:',   '%14.5f%14.5f%14.5f',[1,3], 'box_dcos1'    , []
     'selection dcos2:',   '%14.5f%14.5f%14.5f',[1,3], 'box_dcos2'    , []
     'selection dcos3:',   '%14.5f%14.5f%14.5f',[1,3] , 'box_dcos3'    , []
    };

v6p1_field_defs = cat(1, v6_field_defs, ...
    {
    'reordered toplc(lps, mm): ','%14.5f%14.5f%14.5f',[1,3],'reordered_toplc', []
    'reordered center(lps, mm):','%14.5f%14.5f%14.5f',[1,3],'reordered_center',[]
    'reordered spacing(mm):    ','%14.5f%14.5f%14.5f',[1,3],'reordered_spacing'  , []
    'reordered number of points:',   '%d %d %d',[1,3], 'reordered_n_points'  , []
    'reordered dcos1:',     '%14.5f%14.5f%14.5f',[1,3], 'reordered_dcos1'    , []
    'reordered dcos2:',     '%14.5f%14.5f%14.5f',[1,3], 'reordered_dcos2'   , []
    'reordered dcos3:',     '%14.5f%14.5f%14.5f',[1,3], 'reordered_dcos3'   , []
    });

separator = '===================================================';

% Open file for writing

ddf.filename = strcat(outputname, '.ddf');
outfile = fopen(ddf.filename,'w');

if (ddf.version == 6.1) 
  ddf.version = '6.1';
  field_defs  = v6p1_field_defs;
end

if (ddf.version == 6) 
  ddf.version = '6';
  field_defs  = v6_field_defs;
end

if (ddf.version == 5.1) 
  ddf.version = '5.1';
end

if (ddf.version == 5) 
  ddf.version = '5';
end

% =======================================================================
% Begin writing the actual data (v6)
% =======================================================================

if (strcmp(ddf.version, '6.1') | strcmp(ddf.version, '6'))
  
  number_of_fields = size(field_defs,1);
  
  fprintf(outfile, 'DATA DESCRIPTOR FILE\n');  
  
  for (fieldnum = 1:22)
    if (~isfield(ddf, field_defs{fieldnum,4}))
      error('Unable to find field "%s" in your ddf structure "%s"', ...
            field_defs{fieldnum,4}, inputname(2));
    end
    fprintf(outfile, sprintf('%s %s', ...
                             field_defs{fieldnum, 1}, ...
                             field_defs{fieldnum, 2}), ...
            getfield(ddf, field_defs{fieldnum,4}));
    fprintf(outfile, '\n');
  end
  
  % Information about each dimension
  
  fprintf(outfile, sprintf('dimension 1: type: %s npoints: %d\n', ...
                           ddf.dimension_type{1}, ...
                           ddf.specpoints));
  for (i = 1:3)
    fprintf(outfile, ...
            sprintf(['dimension %d: type: %s ' ...
                    'npoints: %d pixel spacing(mm): %.6f'],...
                    i+1, ...
                    ddf.dimension_type{i+1}, ...
                    ddf.npix(i), ...
                    ddf.pixel_spacing(i)));    
    fprintf(outfile, '\n');
  end
  
  for (fieldnum = 34:38)
    if (~isfield(ddf, field_defs{fieldnum,4}))
      error('Unable to find field "%s" in your ddf structure "%s"', ...
            field_defs{fieldnum,4}, inputname(2));
    end
    fprintf(outfile, sprintf('%s %s', ...
                             field_defs{fieldnum, 1}, ...
                             field_defs{fieldnum, 2}), ...
            getfield(ddf, field_defs{fieldnum,4}));
    fprintf(outfile, '\n');
  end
  
  % BEGIN OUTPUT OF MR PARAMETERS
  
  fprintf(outfile, '%s\n', separator);
  fprintf(outfile, 'MR Parameters\n');

  for (fieldnum = 39:49)
    if (~isfield(ddf, field_defs{fieldnum,4}))
      error('Unable to find field "%s" in your ddf structure "%s"', ...
            field_defs{fieldnum,4}, inputname(2));
    end
    fprintf(outfile, sprintf('%s %s', ...
                             field_defs{fieldnum, 1}, ...
                             field_defs{fieldnum, 2}), ...
            getfield(ddf, field_defs{fieldnum,4}));
    fprintf(outfile, '\n');
  end 
  
  % Sat band information
  
  for (satnum = 1:ddf.num_sat_bands)
    fprintf(outfile, ...
            'sat band %d thickness(mm): %.6f', ...
            satnum, ddf.satband_thickness(satnum));
    fprintf(outfile, '\n');
    fprintf(outfile, ...
            'sat band %d orientation:        %14.5f%14.5f%14.5f', ...
            satnum, ddf.satband_orientation(satnum,:));  
    fprintf(outfile, '\n');
    fprintf(outfile, ...
            'sat band %d position(lps, mm):  %14.5f%14.5f%14.5f', ...
            satnum, ddf.satband_position(satnum,:));  
    fprintf(outfile, '\n');
  end

  % BEGIN OUTPUT OF SPECTROSCOPY PARAMETERS
  
  fprintf(outfile, '%s\n', separator);
  fprintf(outfile, 'Spectroscopy Parameters\n');  

  for (fieldnum = 50:number_of_fields)
    if (~isfield(ddf, field_defs{fieldnum,4}))
      error('Unable to find field "%s" in your ddf structure "%s"', ...
            field_defs{fieldnum,4}, inputname(2));
    end
    fprintf(outfile, sprintf('%s %s', ...
                             field_defs{fieldnum, 1}, ...
                             field_defs{fieldnum, 2}), ...
            getfield(ddf, field_defs{fieldnum,4}));
    fprintf(outfile, '\n');
  end 

  fprintf(outfile, '%s\n', separator);  
  
  % BEGIN OUTPUT OF ANY OTHER PROCESSING RESULTS
  
  if (isfield(ddf, 'misc_info'))
    fprintf(outfile, '%s', ddf.misc_info);
  end
  
end

% =======================================================================
% Begin writing the actual data (v5); this part is different from
% v5 allows multiple data fields to be listed on the same line
% =======================================================================

if (strcmp(ddf.version, '5.1') | strcmp(ddf.version, '5'))

  fprintf(outfile,'------------------------------\n');
  fprintf(outfile,'DATA DESCRIPTOR FILE version 5\n');
  fprintf(outfile,'------------------------------\n');
  fprintf(outfile,'\n');
  fprintf(outfile,'GENERAL PARAMETERS\n');
  fprintf(outfile,'study id: %s\n', ddf.studyid);
  fprintf(outfile,'study # : %s\n', ddf.studynum);
  fprintf(outfile,'position: %s\n', ddf.position);
  fprintf(outfile,'coil: %s\n', ddf.coil);
  fprintf(outfile,'series #: %d\n', ddf.series);
  fprintf(outfile,'orientation: %5d', ddf.orientation);
  if (ddf.orientation == 13)
    fprintf(outfile, '     axial RLxAPxSI');
  end
  fprintf(outfile,'\n');
  fprintf(outfile,'echo/time/met index: %5d     value: %10.2f\n', ...
          ddf.index, ddf.value);
  fprintf(outfile,'root name: %s\n', outputname);
  fprintf(outfile,'comment: %s\n', ddf.comment);  

  fprintf(outfile,'\n');
  fprintf(outfile,'SEQUENCE INFORMATION\n');
  fprintf(outfile,'sequence name: %s\n', ddf.sequence);
  fprintf(outfile,'localization type: %3d\n', ddf.localizationtype);

  if (isfield(ddf, 'satpulsepos'))
    fprintf(outfile,'spec frequency(MHz): %10.2f\n', ddf.specfrequency);
    fprintf(outfile,'sweepwidth(Hz):  %10.2f    dwelltime: %10.1f\n', ...
            ddf.sweepwidth, ddf.dwelltime);
    fprintf(outfile,'satpulse pos(Hz):  %10.2f    bandwidth(Hz): %10.1f\n', ...
            ddf.satpulsepos, ddf.bandwidth);  
  else
    fprintf(outfile,'spec frequency(kHz): %10.5f \n', ...
            ddf.specfrequency);
    fprintf(outfile,'sweepwidth(Hz):  %9.1f     dwelltime: %10.1f\n', ...
            ddf.sweepwidth, ddf.dwelltime);
    fprintf(outfile,'centfreq pos(Hz): %10.1f       pulse type: %8d\n', ...
            ddf.centfreq, ddf.pulsetype);
  end
    
  fprintf(outfile,'beginning of acquisition (ms): %10.3f\n', ddf.acqbeginning);
  
  fprintf(outfile, 'gradient on time (ms): %10.3f\n', ...
        ddf.gradontime);
  fprintf(outfile, ...
          'nex: %6d     chop: %3d     even_sym: %3d     pe_order: %3d\n', ...
          ddf.nex, ddf.chop, ddf.even_sym, ddf.pe_order);
  
  fprintf(outfile,'\n');
  fprintf(outfile,'DATA FORMAT INFORMATION\n');
  fprintf(outfile,'datatype: %5d     floating point complex format\n', ...
          ddf.datatype);
  fprintf(outfile,'number of dimensions: %3d\n', ...
          ddf.numdim);
  
  fprintf(outfile,'DIM:  1  type 1  npts: %5d  sw(Hz): %8.1f  acqpts: %5d\n', ...
          ddf.specpoints, ddf.sweepwidth, ddf.acqspec);
  fprintf(outfile,'DIM:  2  type 3  npts: %5d    fov:  %8.1f  acqpts: %5d\n', ...
          ddf.npix(1), ddf.fov(1), ddf.acqpts(1));
  fprintf(outfile,'DIM:  3  type 3  npts: %5d    fov:  %8.1f  acqpts: %5d\n', ...
          ddf.npix(2), ddf.fov(2), ddf.acqpts(2));
  fprintf(outfile,'DIM:  4  type 3  npts: %5d    fov:  %8.1f  acqpts: %5d\n', ...
          ddf.npix(3), ddf.fov(3), ddf.acqpts(3));
  
  fprintf(outfile,'\n');
  fprintf(outfile,'RECONSTRUCTED PARAMETERS IN LPS COORDINATES\n');
  fprintf(outfile,'center: %14.5f%14.5f%14.5f\n', ddf.recon_center);
  fprintf(outfile,'toplc:  %14.5f%14.5f%14.5f\n', ddf.recon_toplc);
  fprintf(outfile,'dcos:   %14.5f%14.5f%14.5f\n', ddf.recon_dcos1);
  fprintf(outfile,'dcos:   %14.5f%14.5f%14.5f\n', ddf.recon_dcos2);
  fprintf(outfile,'dcos:   %14.5f%14.5f%14.5f\n', ddf.recon_dcos3);
  fprintf(outfile,'PHASE ENCODE PARAMETERS IN LPS COORDINATES\n');
  fprintf(outfile,'pcenter:%14.5f%14.5f%14.5f\n', ddf.encode_center);
  fprintf(outfile,'pfov:   %14.5f%14.5f%14.5f\n', ddf.encode_fov);
  fprintf(outfile,'pmatrix:%14d%14d%14d\n', ddf.encode_matrix);
  fprintf(outfile,'reverse:%14d%14d%14d\n', ddf.encode_reverse);
  fprintf(outfile,'pdcos:  %14.5f%14.5f%14.5f\n', ddf.encode_dcos1);
  fprintf(outfile,'pdcos:  %14.5f%14.5f%14.5f\n', ddf.encode_dcos2);
  fprintf(outfile,'pdcos:  %14.5f%14.5f%14.5f\n', ddf.encode_dcos3);
  fprintf(outfile,'BOX PARAMETERS IN LPS COORDINATES\n');
  fprintf(outfile,'bcenter:%14.5f%14.5f%14.5f\n', ddf.box_center);
  fprintf(outfile,'bsize:  %14.5f%14.5f%14.5f\n', ddf.box_size);
  fprintf(outfile,'bdcos:  %14.5f%14.5f%14.5f\n', ddf.box_dcos1);
  fprintf(outfile,'bdcos:  %14.5f%14.5f%14.5f\n', ddf.box_dcos2);
  fprintf(outfile,'bdcos:  %14.5f%14.5f%14.5f\n', ddf.box_dcos3);

  fprintf(outfile,'refacq:   %s\n',ddf.refacq);
  fprintf(outfile,'refrecon: %s\n',ddf.refrecon);
  
end

fclose(outfile);

