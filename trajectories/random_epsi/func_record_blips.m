function record = func_record_blips(loc_samp_3d, type_pattern, display)

%
% blip_record = func_record_blips_v2(loc_samp_3d, encodeMatrix, 1);
%
% Creates a structure that has x encode by y encode entries. Associated
% with each entry is a value indicating type of sampling pattern (copy of
% encodeMatrix) and associated x/y blips with appropriate amplitude. Also the 
% total number of TRs is a field. There is also an option to display these blips.
%
% Inputs: a 3d sampling matrix created by func_create3d_v2.m, the 2D encode
% matrix with the type of blips for each k-space patch, and a variable
% indicating if the blip patterns should be displayed (0 - not displayed, 
% 1 - displayed).
%
% Output: resulting structure identifiying type of pattern and blip
% sequence for each phase encode step and numTRs.

% Simon Hu 11/26/07
% modified v2 9/1/08

[length_x length_y length_f] = size(loc_samp_3d);

% create matrices to track type and blip patterns
xblip_pattern = zeros(length_x,length_y,length_f);
yblip_pattern = zeros(length_x,length_y,length_f);

% step through the phase x by y phase encode grid and determine blips for
% each. This uses the convention in func_create3d_v2.m
for x = 1:length_x
    for y = 1:length_y
        switch (type_pattern(x,y))
            case 2
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,2,1,length_f);
            case 3
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,1,2,length_f);
            case 4
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,2,2,length_f);
            case 5
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,3,1,length_f);
            case 6
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,1,3,length_f);
            case 7
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,3,2,length_f);
            case 8
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,2,3,length_f);
            case 9
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,3,3,length_f);
            case 10
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,4,1,length_f);
            case 11
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,1,4,length_f);
            case 12
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,4,2,length_f);
            case 13
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,4,3,length_f);                
            case 14
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,2,4,length_f);
            case 15
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,3,4,length_f);
            case 16
                [xblip_pattern yblip_pattern] = fill_with_blips(xblip_pattern,yblip_pattern,loc_samp_3d,x,y,4,4,length_f);
        end
    end
end

% figure out number of TRs 
numTRs = sum(type_pattern(:) > 0);

% display blips
if (display == 1)
    figure;  
    for y = 1:length_y
        for x = 1:length_x
            if (type_pattern(x,y) > 0)
                subplot(2,1,1);
                plot_blip_pattern(xblip_pattern,x,y);
                title(sprintf('viewx = %d, viewy = %d - xblip pattern',x,y));
                xlabel('x blip before nth flyback plateau');
                subplot(2,1,2);
                plot_blip_pattern(yblip_pattern,x,y);
                title('yblip pattern');
                xlabel('y blip before nth flyback plateau');
                pause;
            end
        end
    end 
end
                
% put generated data into a struct
record = struct('type_of_blip_pattern', type_pattern, 'x_blip_pattern', xblip_pattern, 'y_blip_pattern', yblip_pattern, 'numTRs', numTRs);

% --------------------Subfunctions-----------------------------------------

function [xp yp] = fill_with_blips(xp_in,yp_in,matrix3d,x,y,xpsan,yspan,length_f) 
% Fills in the x/y blip patterns associated with a patch of k-space.
% Determines the size of blip needed to move to next k-space point.
xp = xp_in;
yp = yp_in;
% Assign blips
prev_xoff = 1;
prev_yoff = 1;
for f = 1:length_f
    [curr_xoff curr_yoff] = find(squeeze(matrix3d(x:x+xpsan-1,y:y+yspan-1,f)) == 1);
    xp(x,y,f) = prev_xoff - curr_xoff;
    yp(x,y,f) = prev_yoff - curr_yoff;
    prev_xoff = curr_xoff;
    prev_yoff = curr_yoff;
end

function plot_blip_pattern(blip_pattern,x,y)
% Plots triangular blips of -1,0,1 polarity according to input sequence of
% numbers.
len = length(blip_pattern(x,y,:));
N = 50;
pattern = zeros(1,2*len*N);

for f = 1:len
    pattern((f-1)*2*N+1:(f-1)*2*N+N) = blip_pattern(x,y,f)*triangle(1,N);
end
blip_n = ((0:length(pattern)-1) + 3*N/2)/(2*N);
plot(blip_n,pattern);
ylim([-3 3]);

function t = triangle(G,N)
% Returns a triangle of length N and height G. 
if (rem(N,2) == 1)
    mid_N = round(N/2);
    t(1:mid_N) = 0:G/(mid_N-1):G;
    t(mid_N:N) = G:-G/(mid_N-1):0;
else
    mid_N = N/2;
    t(1:mid_N) = 0:G/(mid_N-1):G;
    t(mid_N+1:N) = G:-G/(mid_N-1):0;
end
