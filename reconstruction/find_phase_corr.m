function phase = find_phase_corr(S)
% phase = find_phase_corr(S)
% simple way that may work to find zero-order phase correction
%
% estimates phase to put the most signal in the real channel

% weight sum by magnitude
phase = -angle( sum(S(:) .* abs(S(:))) );