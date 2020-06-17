
% rawSpectra must be an N x M matrix (N = frequency points, M = receive array elements)
% noiseMask must be an N x 1 logical matrix (N = frequency points).
%           true --> this point is baseline noise.
%           false --> this point contains signal.

function [svdCoilAmplitudes,svdWeights,svdRecombination, svdQuality] = svdRecon(rawSpectra, noiseMask, varargin)

if nargin<3
    options = struct();
elseif nargin==3
    options = varargin{1};
else
    % Wrap any cell arrays before passing in to struct()
    for idx=2:2:numel(varargin)
        if iscell(varargin{idx})
            varargin{idx} = {varargin{idx}};
        end
    end
    
    options = struct(varargin{:});
end

if ~isfield(options,'debug')
    options.debug = 0;
end

if ~isfield(options,'phaseRefChannel')
    options.phaseRefChannel = 1;
end

if numel(size(rawSpectra))>2
    error('rawSpectra must be an N x M matrix (N = frequency points, M = receive array elements)')
end

%% Allocate storage for some results
%sizeSpectra = size(rawSpectra,1);
nCoils = size(rawSpectra,2);

%% Measure noise statistics from data (unless instructed otherwise)
if isfield(options,'noiseCov') && ~ischar(options.noiseCov)
  % Noise covariance matrix has been supplied
  noiseCov=options.noiseCov;
elseif isfield(options,'noiseCov') && strcmp(options.noiseCov,'disable')
    % Disable noise prewhitening
    noiseCov=eye(nCoils)*0.5;
elseif isfield(options,'noiseCov') && strcmp(options.noiseCov,'diag')
    % Use only the noise variances (not off diagonal elements)...
    noiseCov=diag(diag(cov(rawSpectra(noiseMask,:))));
elseif ~isfield(options,'noiseCov') || strcmp(options.noiseCov,'normal')
    % Estimate from data (DEFAULT)
    noiseCov=cov(rawSpectra(noiseMask,:));
else
    error('Option "noiseCov" has an unknown value.')
end

[noiseVec, noiseVal] = eig(noiseCov);

% (This formula is a simple adaptation of that in
% testNoiseUncorrelationNmrFft2.m, adding sim.imagingFrequency to the mix.)
scaleMatrixFft = noiseVec*diag(sqrt(0.5)./sqrt(diag(noiseVal)));
invScaleMatrixFft = inv(scaleMatrixFft);

scaledSpectra=rawSpectra*scaleMatrixFft;

if options.debug
    disp('DEBUG: Please confirm that this matrix is approximately 0.5*eye(nCoils)')
    cov(scaledSpectra(noiseMask,:))
    
    disp('Noise covariance matrix eigenvalues:')
    disp(noiseVal)

    figure(9);clf;mypcolor(1:nCoils,1:nCoils,abs(noiseCov));colorbar
    set(gca,'YTick',1:nCoils,'YTickLabel',options.coilNames,...
        'XTick',1:nCoils,'XTickLabel',options.elementNames)
end

%% Compute optimal reconstruction using SVD
[u,s,v]=svd(scaledSpectra,'econ');
% SVD quality indicator
svdQuality(1) = ((s(1,1)/norm(diag(s)))*sqrt(nCoils)-1)/(sqrt(nCoils)-1);

% Coil amplitudes
svdCoilAmplitudes=v(:,1)'*invScaleMatrixFft;
% There's an arbitrary scaling here such that the first coil weight is
% real and positive

svdRescale = norm(svdCoilAmplitudes)*normalise(svdCoilAmplitudes(options.phaseRefChannel));
%svdRescale = norm(svdCoilAmplitudes);

svdCoilAmplitudes=svdCoilAmplitudes / svdRescale;

if options.debug
    fprintf('DEBUG: SVD quality indicator = %g\n',svdQuality);
    fprintf('SVD coil amplitudes:\n');
    disp(svdCoilAmplitudes)
end

svdRecombination = u(:,1)*s(1,1)*svdRescale;

%% Added 2 Feb 2014:
% svdWeights = 0.5*svdCoilAmplitudes'*inv(noiseCov) * conj(svdRephase) * svdRephase; % V4 version. Fails. Wrong matrix sizes.

svdWeights = 0.5*inv(noiseCov) * svdCoilAmplitudes' * conj(svdRescale) * svdRescale; % Try this. Not derived carefully. But tested numerically and it gives the same answer to within 1e-15.
svdWeights = permute(svdWeights,[2 1]);         
%             nawmfilename = sprintf('%s/%s/%s/spectra/%s_wm_pmask',path_dir,pathname_spectra_sw{nptn},bnum,tnum,tnum);
%             t2allfilename = sprintf('%s/%s/%s/spectra/%s_t2all_pmask',path_dir,pathname_spectra_sw{nptn},bnum,tnum,tnum);

