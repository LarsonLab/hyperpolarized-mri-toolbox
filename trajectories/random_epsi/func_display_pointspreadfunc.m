function psf_ratio = func_display_pointspreadfunc(loc_samp_3d,type)

%
% ratio = func_diplay_pointspreadfunc(loc_samp_3d,4);
%
% Shows a graphical representation of the point spread function associated 
% with a 3d trajectory. If type 1,2,3 are selected successive slices along 
% x,y,and f respectively are displayed instead. If type 4 is selected, a 3d
% scatter plot is displayed. Also the peak to maximum sidelobe ratio of the
% psf is output for every case. If none of the valid types 1-4 is selected,
% e.g. selected 0 or 5, then nothing is displayed and only the peak to max
% sidelobe ratio is returned.
%
% Input: 3d trajectory matrix, type of display
%
% Output: peak-to-max-sidelobe ratio of point spread function

% Simon Hu 11/25/07

[length_x length_y length_f] = size(loc_samp_3d);

% take 3d fourier transform to get psf

psf = fftshift(fft(loc_samp_3d,[],1),1);
psf = fftshift(fft(psf,[],2),2);
psf = fftshift(fft(psf,[],3),3);
psf = abs(psf);
psf = psf/max(max(max(psf)));

% display psf - type 1,2,3 are multiple 2d plots, 4 is a 3d scatter plot

if (type == 1) 
    
    figure;
    for x = 1:length_x
        image(squeeze(255*psf(x,:,:)));
        title(sprintf('x = %d',x)); xlabel('f'); ylabel('y');
        pause;
    end
    
elseif (type == 2)
    
    figure;
    for y = 1:length_y
        image(squeeze(255*psf(:,y,:)));
        title(sprintf('y = %d',y)); xlabel('f'); ylabel('x');
        pause;
    end
    
elseif (type == 3)
    
    figure;
    for f = 1:length_f
        image(squeeze(255*psf(:,:,f)));
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
    % take a zero value, so eps needs to be added
    scatter3(x,y,f,100*psf(:)+eps,'filled');
    title('3d point spread function'); xlabel('x'); ylabel('y'); zlabel('f');
    
else
    
    % no operation
    temp = loc_samp_3d;
    
end

% report psf peak to max sidelobe ratio

psf_list = psf(:);
ind = find(psf_list == max(psf_list));
peak = psf_list(ind);
psf_list(ind) = 0;
ind2 = find(psf_list == max(psf_list));
max_sidelobe = psf_list(ind2(1));
psf_ratio = peak/max_sidelobe;
