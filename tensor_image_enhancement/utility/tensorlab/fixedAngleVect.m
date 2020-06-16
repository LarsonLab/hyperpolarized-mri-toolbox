function v = fixedAngleVect(N, d, angle, randomize)
%FIXEDANGLEVECT Vectors with fixed angle.
%    V = FIXEDANGLEVECT(N, d, angle) computes d vectors of
%    length N. Each vector in v has unit norm and a fixed angle with all
%    other vectors. d must be lower than or equal to N + 1,
%    otherwise, only N + 1 vectors are returned. The angle should be given
%    in radians.
%
%    V = FIXEDANGLEVECT(N, d, angle, randomize) multiplies the resulting
%    vectors with a random unitary matrix if randomize is true. The default
%    is true. 

% References:
% [1] Simplex - http://en.wikipedia.org/wiki/Simplex (visited 2013/04/17)
% 
% Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%          Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2014/06/10   NV      Fixed bug that the cosine of the angle instead of
%                        the angle is required. 
% - 2013/04/17   NV      Initial version


    if nargin < 4, randomize = true; end
    
    % checks
    d = min(N + 1, d);

    % initialization
    v = zeros(N, d);

    % Algorithm
    for k = 1:d - 1
        v(k, k) = sqrt(1 - sum(v(1:k-1,k).^2));
        v(k, k+1:end) = (cos(angle) - sum(v(1:k-1,k).^2))./ v(k, k);
    end
    v(d, d) = sqrt(1 - sum(v(1:d-1,d).^2));

    % Randomize
    if randomize
        rot = orth(randn(N));
        v = rot*v;
    end
end