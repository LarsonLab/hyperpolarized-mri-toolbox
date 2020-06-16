function C4 = stcum4(Y, L)
%STCUM4 Fourth-order spatio-temporal cumulant tensor.
%   C4 = stcum4(Y,L) computes fourth-order spatio-temporal cumulant c4 of 
%   the input matrix Y in which each row is an observation and each column 
%   is a variable. Herein,
%
%      C4(i,j,k,l,t1,t2,t3) = E[yi(t).*conj(yj(t+t1)).*conj(yk(t+t2)).*yl(t+t3)] ...
%                    - E[yi(t).*conj(yj(t+t1))]*E[conj(yk(t+t2)).*yl(t+t3)] ...
%                    - E[yi(t).*conj(yk(t+t2))]*E[conj(yj(t+t1)).*yl(t+t3)] ...
%                    - E[yi(t).*yl(t+t3)]*E[conj(yj(t+t1)).*conj(yk(t+t2))]
%
%   where the expectation E is approximated by the arithmetic mean and xi
%   is the i-th mean centered variable, Y(:,i)-mean(Y(:,i)) (and
%   analogously for xj, xk and xl). The variables t1, t2 and t3 can vary
%   between -L:L and are indexed as 1:2L+1.
%
%   C4 = stcum4(Y) is equivalent to stcum4(Y,0).
%
%   See also xcum4.

%   Authors: Frederik Van Eeghem (Frederik.VanEeghem@esat.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

if nargin < 2
    L = 0;
end

M = size(Y,2);

% Center data
Y = bsxfun(@minus, Y, mean(Y));

% Compute spatio-temporal cumulant by computing the cross cumulant for all
% time shifts.
C4 = zeros([M M M M (2*L+1) (2*L+1) (2*L+1)]);
t = 0:size(Y,1)-2*L-1;
if isreal(Y)
    for t1 = 1:2*L+1
        for t2 = t1:2*L+1
            for t3 = t2:2*L+1
                tmp = xcum4(Y(t+L+1,:),Y(t+t1,:),Y(t+t2,:),Y(t+t3,:),0);
                C4(:,:,:,:,t1,t2,t3) = tmp;
                C4(:,:,:,:,t2,t1,t3) = permute(tmp,[1 3 2 4]);
                C4(:,:,:,:,t2,t3,t1) = permute(tmp,[1 3 4 2]);
                C4(:,:,:,:,t1,t3,t2) = permute(tmp,[1 2 4 3]);
                C4(:,:,:,:,t3,t1,t2) = permute(tmp,[1 4 2 3]);
                C4(:,:,:,:,t3,t2,t1) = permute(tmp,[1 4 3 2]);
            end
        end
    end
else 
    for t1 = 1:2*L+1
        for t2 = t1:2*L+1
            for t3 = 1:2*L+1
                tmp = xcum4(Y(t+L+1,:),Y(t+t1,:),Y(t+t2,:),Y(t+t3,:),0);
                C4(:,:,:,:,t1,t2,t3) = tmp;
                C4(:,:,:,:,t2,t1,t3) = permute(tmp,[1 3 2 4]);
            end
        end
    end
end

end