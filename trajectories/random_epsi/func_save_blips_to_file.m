function func_save_blips_to_file(filename, blip_record)

%
% func_save_blips_to_file('blip_file', blip_record);
% 
% Saves blip patterns to a text file readable in C. The first four numbers
% are the x resolution, y resolution, f resolution, and actual number of TRs
% to be played out. Then x*y*f values for the xblips are written, followed 
% by the x*y*f yblip values. Lastly, blip pattern type for each view is 
% written, letting EPIC know which views  to skip (0 is skip). All values 
% are separated by a space.
%
% Inputs: file name, blip record created by func_record_blips.m
%
% Output: Creates a text file where the blip record is written. 

% Simon Hu 2/18/08

fid = fopen(filename, 'wt');
xblips = blip_record.x_blip_pattern;
yblips = blip_record.y_blip_pattern;
type_blips = blip_record.type_of_blip_pattern;

% Write dimensions and numTRs to file.

[length_x length_y length_f] = size(xblips);
numTRs = blip_record.numTRs;
fprintf(fid, '%d %d %d %d ', length_x, length_y, length_f, numTRs);

% Write x*y*f xblip values. The order is the f points for viewx=1,viewy=1, 
% then the f points for viewx=2,viewy=1, and so forth. 

for y = 1:length_y
    for x = 1:length_x
        fprintf(fid, '%d ', xblips(x,y,:));
    end
end

% Write x*y*f yblip values. Same write order.

for y = 1:length_y
    for x = 1:length_x
        fprintf(fid, '%d ', yblips(x,y,:));
    end
end

% Write type values to let EPIC know which views to skip. Same write order.

for y = 1:length_y
    for x = 1:length_x
        fprintf(fid, '%d ', type_blips(x,y));
    end
end

fclose(fid);
