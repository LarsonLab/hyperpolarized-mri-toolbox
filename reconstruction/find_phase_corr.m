function phase = find_phase_corr(S)
% simple way that may work to find zero-order phase correction
%

% weight sum by magnitude
phase = -angle( sum(S(:) .* abs(S(:))) );