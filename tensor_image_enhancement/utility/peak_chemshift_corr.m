% automated peak finding algorithm
% v3. based on cross correlation between voxel of interest and a standard
% voxel


function [notfound closest_Idx] = peak_chemshift_corr(a1_standard, a1_signal, empirical_center, varargin)
% a1 is the interpolated spectrum
% a1_standard is the voxel used as example
% a1_signal is the voxel of interest

varargs_length = 2;
if length(varargin) > varargs_length
    error('myfuns:peak_chemshift:TooManyInputs', ...
        'requires at most %d optional inputs',varargs_length);
end

% optargs = {detection_method, detection_direction, threshold_peakshift,
% AUC_width}
optargs = {[500],[100]}; % defaults
for i_args = 1:length(varargin),
    if ~isempty(varargin{i_args})
        optargs{i_args} = varargin{i_args};
    end
end

notfound = false;

% the width to which cross correlation was calculated
% roughly equal to pyruvate peak width
options.width_search = optargs{1};
options.chemshift_range = optargs{2};

% ---- calculate cross correlation ----
% between standard spectrum vs chemshifted spectrum
a1_standard_seg = a1_standard(empirical_center-options.width_search/2:empirical_center+options.width_search/2);
a1_signal_seg = a1_signal(empirical_center-options.width_search/2:empirical_center+options.width_search/2);

a1_corr = xcorr(a1_signal_seg,a1_standard_seg);
[max_temp, I] = max(a1_corr,[],2);
closest_Idx = I-1-options.width_search + empirical_center;

% sanity check, if the peak is way off, then probably wrong
% use empirical center by default in this case
if abs(closest_Idx-empirical_center) > options.chemshift_range/2;
    notfound = true;
    closest_Idx = empirical_center;
end

    
end