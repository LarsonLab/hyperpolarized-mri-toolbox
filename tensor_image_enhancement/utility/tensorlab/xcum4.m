function [c4,m4] = xcum4(X1,X2,X3,X4,center)
%XCUM4 Fourth-order cross-cumulant tensor.
%   [c4,m4] = xcum4(X1,X2,X3,X4) computes fourth-order cross moment m4 and
%   fourth-order cross cumulant c4 of 4 matrices X1 to X4 in which each row
%   is an observation and each column is a variable. Herein,
%
%      m4(i,j,k,l) = E[x1i.*conj(x2j).*conj(x3k).*x4l]
%      c4(i,j,k,l) = E[x1i.*conj(x2j).*conj(x3k).*x4l] ...
%                    - E[x1i.*conj(x2j)]*E[conj(x3k).*x4l] ...
%                    - E[x1i.*conj(x3k)]*E[conj(x2j).*x4l] ...
%                    - E[x1i.*x4l]*E[conj(x2j).*conj(x3k)]
%
%   where the expectation E is approximated by the arithmetic mean and xri
%   is the i-th mean centered variable, Xr(:,i)-mean(Xr(:,i)) (and
%   analogously for xrj, xrk and xrl).
%
%   [c4,m4] = xcum4(X1,X2,X3,X4,center) controls the centering of the 
%   input columns so that their mean is zero. The default value for
%   the center argument is true, meaning the columns will be centered.
%   Centering can be turned off by setting the center argument to false or
%   'nocenter'.
%
%   See also cum4.

%   Authors: Frederik Van Eeghem (Frederik.VanEeghem@esat.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the center option.
if nargin < 5, center = true; end
if ischar(center), center = ~strcmpi(center,'nocenter'); end

[N1,M1] = size(X1);
[N2,M2] = size(X2);
[N3,M3] = size(X3);
[N4,M4] = size(X4);

Nvec = [N1,N2,N3,N4];
if any(Nvec ~= N1)
    error('xcum4:inputs','Inputs must all have the same number of rows.');
end

% Center variables
if center
    X1 = bsxfun(@minus, X1, mean(X1));
    X2 = bsxfun(@minus, X2, mean(X2));
    X3 = bsxfun(@minus, X3, mean(X3));
    X4 = bsxfun(@minus, X4, mean(X4));
end

% Compute the fourth-order cross moment
% m4(i,j,k,l) = E[x1i*conj(x2j)*conj(x3k)*x4l].
temp1 = kr(X2',X1.');
temp2 = kr(X4.',X3');
m4 = temp1*temp2.'/N1;

% Compute second order contributions.
nv = ones(N1,1)/N1;
Y1 = temp1*nv;
Y2 = temp2*nv;
Y3 = X1.'*conj(X3)/N1;
Y4 = X2'*X4/N1;
Y5 = X1.'*X4/N1;
Y6 = X2'*conj(X3)/N1;
temp = kron(Y6, Y5);
temp = reshape(temp,[M1*M2,M4,M3]);
temp = permute(temp,[1,3,2]);
temp = reshape(temp,[M1*M2,M3*M4]);

% Combine the information to obtain the matricized cumulant.
c4 = m4 - Y1*Y2.' - kron(Y4, Y3) - temp;

% Tensorize the moment and cumulant.
c4 = reshape(c4,[M1 M2 M3 M4]);
m4 = reshape(m4,[M1 M2 M3 M4]);

end