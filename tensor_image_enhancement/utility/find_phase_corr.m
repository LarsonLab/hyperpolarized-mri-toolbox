function phase = find_phase_corr(S,varargin)
% v0. by Peder Larson
% v1. add multidimensional feature (timepoints, metabolites etc.) 20170821
% S(f,voxel,time,etc)

varargs_length = 1;
if length(varargin) > varargs_length
    error('myfuns:find_phase_corr:TooManyInputs', ...
        'requires at most %d optional inputs',varargs_length);
end

optargs = {[]}; % defaults
for i_args = 1:length(varargin),
    if ~isempty(varargin{i_args})
        optargs{i_args} = varargin{i_args};
    end
end

options.SNR_threshold = optargs{1};

% weight sum by square of magnitude
phase = -angle( squeeze(sum(S .* abs(S),1)) );


end