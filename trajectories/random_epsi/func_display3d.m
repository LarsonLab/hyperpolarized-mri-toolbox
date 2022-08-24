function func_display3d(loc_samp_3d,type)

%
% func_display3d(loc_samp_3d,4);
%
% Displays slice by slice the undersampling pattern associated with a 3d
% trajectory or the entire sampling pattern as a 3d scatter plot.
%
% Inputs: 3d matrix of zeros and ones and a type number that specifies
% which axis (x = 1, y = 2, f = 3, 3d = 4) over which to step through. 
% For example, selecting type = 2 will display x-f slices associated with 
% every y point. Selecting type 4 will produce a single 3d scatter plot.

% Simon Hu 11/24/07

[length_x length_y length_f] = size(loc_samp_3d);

if (type == 1) 
    
    figure;
    for x = 1:length_x
        imagesc(squeeze(loc_samp_3d(x,:,:)));
        title(sprintf('x = %d',x)); xlabel('f'); ylabel('y');
        pause;
    end
    
elseif (type == 2)
    
    figure;
    for y = 1:length_y
        imagesc(squeeze(loc_samp_3d(:,y,:)));
        title(sprintf('y = %d',y)); xlabel('f'); ylabel('x');
        pause;
    end
    
elseif (type == 3)
    
    figure;
    for f = 1:length_f
        imagesc(squeeze(loc_samp_3d(:,:,f)));
        title(sprintf('f = %d',f)); xlabel('y'); ylabel('x');   
        pause;
    end
    
elseif (type == 4)
    
    % 3d scatter plot: create vectors for x,y,z coordinates for every
    % single point in the 3d matrix. Basically have to specify value for
    % every coordinate to use the scatter3 function.
    
    figure;
    len = length(loc_samp_3d(:));
    
    % for x, repeat x vector length_y*length_f times
   
    x = zeros(1,len);
    for a = 1:length_x
        x(a:length_x:len) = a;
    end
            
    % for y, repeat y vector length_f times but each value of y has
    % multiplicity of length_x, e.g. if length x is 2, y would go
    % 1,1,2,2,3,3,...1,1,2,2,3,3,...
    
    y = zeros(1,len);
    for a = 1:length_y
        for b = 1:length_x
            start = (a-1)*length_x + b;
            y(start:(length_x*length_y):len) = a;
        end
    end
    
    % for f, write f vector once but with multiplicity length_x*length_y,
    % e.g. if length x is 2 and length y is 4, f would go
    % 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,....
    
    f = zeros(1,len);
    for a = 1:length_f
        start = (a-1)*length_x*length_y + 1;
        f(start:(start+length_x*length_y-1)) = a;
    end
   
    % scatter3: takes x,y,z vectors and plots a dot. note that it cannot
    % take a zero value, so eps needs to be added.
    S = 50000/len; 
    scatter3(x,y,f,S*loc_samp_3d(:)+eps,'filled');
    title('3d sampling pattern'); xlabel('x'); ylabel('y'); zlabel('f');
    
else
    
    % no operation
    temp = loc_samp_3d;
    
end