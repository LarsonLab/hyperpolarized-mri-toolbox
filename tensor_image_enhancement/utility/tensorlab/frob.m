function f = frob(T,squared)
%FROB Frobenius norm
%   FROB(T) returns the Frobenius norm of the tensor T. In the case of
%   incomplete tensors, the weighted Frobenius norm error is used, i.e., the
%   norm of the known elements is computed. In the case T is an efficient
%   representation of a structured tensor, the norm is computed without
%   expanding the tensor.
%
%   FROB(T,'squared') returns the squared Frobenius norm of the tensor T.
%
%   See also norm.

% Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%          Otto Debals         (Otto.Debals@esat.kuleuven.be)
%          Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%          Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%          Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/02/02    OD      Added support Hankel and Loewner
% - 2016/12/18    NV      Added structured tensor support 

if nargin < 2, squared = false; end
switch getstructure(T)
    case 'full'
        f = full(T(:)'*T(:));
    case {'sparse', 'incomplete'}
        f = T.val(:)'*T.val(:);
    case 'cpd'
        W = cellfun(@(u) u'*u, T, 'UniformOutput', false);
        W = prod(cat(3, W{:}),3);
        f = abs(sum(W(:)));
    case 'lmlra'
        W = cellfun(@(u) u'*u, T{1}, 'UniformOutput', false);
        S = tmprod(conj(T{2}), W, 1:length(W), 'T').*T{2};
        f = abs(sum(S(:)));
    case 'btd'
        f = 0;
        for r = 1:length(T)
            W = cellfun(@(u) u'*u, T{r}(1:end-1), 'UniformOutput', false);
            S = tmprod(conj(T{r}{end}), W, 1:length(W), 'T').*T{r}{end};
            f = f + sum(S(:));
            for s = r+1:length(T)
                W = cellfun(@(u,v) u'*v, T{r}(1:end-1), T{s}(1:end-1), ...
                    'UniformOutput', false);
                S = tmprod(conj(T{r}{end}), W, 1:length(W), 'T').*T{s}{end};
                f = f + 2*real(sum(S(:)));
            end
        end
        f = abs(f);
    case 'tt'
        size_tens = cellfun('size', T, 2);
        ttr = [1 cellfun('size', T(2:end), 1) 1];
        tmp = T{1}'*T{1};
        f = tmp(:).';
        for n = 2:length(T)-1
            tmp = reshape(permute(T{n},[2 1 3]),size_tens(n),ttr(n)*ttr(n+1));
            tmp = tmp'*tmp;
            tmp = reshape(tmp, ttr([n n+1 n n+1]));
            tmp = permute(tmp, [1 3 2 4]);
            tmp = reshape(tmp, ttr(n)^2, ttr(n+1)^2);
            f = f * tmp;
        end
        tmp = conj(T{end}*T{end}');
        f = abs(f*tmp(:));
    case 'hankel'
        N = size(T.val,1);
        size_hankel = T.subsize.hankel;
        if T.order==2,
            m = min(T.ind,N-T.ind+1);
            w = [1:m-1 m*ones(1,N-2*m+2) m-1:-1:1];
        else
            w = fft(ones(size_hankel(1),1),N);
            w = w(:);
            for i = 2:numel(size_hankel)
                w = w.*fft(ones(size_hankel(i),1),N);
            end
            w = round(ifft(w).');
        end
        f = w*sum(T.val.*conj(T.val),2);
        
    case 'loewner'
        if T.order > 2
            % full case, as the structured case is not yet supported
            L = ful(T);
            f = frob(L,'squared');
        else
            
            f = T.val(T.ind{1},:);
            g = T.val(T.ind{2},:);
            
            if T.isequidistant
                v = T.structure.v.*conj(T.structure.v);
                tmp = ifft(fft(v,[],1).*fft(sum(f.*conj(f),2),numel(v),1));
                term1 = sum(tmp(end-numel(v)+size(f,1):end,:),1);
                
                tmp = ifft(fft(v(end:-1:1),[],1).*fft(sum(g.*conj(g),2),numel(v),1));
                term2 = sum(tmp(end-numel(v)+size(g,1):end,:),1);
                
                tmp = ifft(bsxfun(@times,fft(v,[],1),fft(f,numel(v),1)),[],1);
                term3 = sum(sum(real(tmp(end-numel(v)+size(f,1):end,:).*conj(g))));
            else
                XY = 1./bsxfun(@minus,T.t(T.ind{1}),T.t(T.ind{2}).');
                XYt2 = (XY.*conj(XY)).';
                term1 = sum(XYt2*sum(f.*conj(f),2));
                term2 = sum(sum(g.*conj(g),2).'*XYt2);
                term3 = sum(sum(real((XYt2*f).*conj(g))));
            end
            
            f = term1+term2-2*term3;
        end
    otherwise
        error('frob:notImplemented', 'Not yet implemented');
end
if ~ischar(squared), f = sqrt(f); end
end