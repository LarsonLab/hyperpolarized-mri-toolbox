function f = frobbtdres(T, U, varargin)
%FROBBTDRES Frobenius norm of residual for a BTD.
%   F = FROBBTDRES(T, U) computes the Frobenius norm of the residual tensor E
%   which is given by btdgen(U)-T.
%
%   FROBBTDRES(T,U,options) or FROBBTDRES(T,U,'key',value) can be used to set
%   the following options:
%
%   - ExpandLimit           Expands structured tensors to full tensors if the
%                           total number of entries prod(getsize(T)) is
%                           smaller than the value of this limit. Default:
%                           1e6.
%   - Type                  The type of the tensor. If not given, the type is
%                           determined automatically. 
%
%   All other options are passed on to BTDRES if needed.
%
%   See also btdres, frob, getstructure
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/30   NV      Initial version
    
    p = inputParser;
    p.addOptional('ExpandLimit', 1e6);
    p.addOptional('Type', getstructure(T));
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    passopt = p.Unmatched;
    
    type = options.Type;
    if any(strcmpi(type, {'full', 'incomplete', 'sparse'}))
        f = frob(btdres(T, U, passopt));
    elseif prod(getsize(T)) <= options.ExpandLimit
        if ~isfield(passopt, 'Format'), passopt.Format = false; end
        f = frob(btdres(ful(T),U, passopt));
    else 
        f = abs(frob(T,'squared') - 2*real(inprod(T,U)) + frob(U,'squared'));
        f = sqrt(f);
    end
end
