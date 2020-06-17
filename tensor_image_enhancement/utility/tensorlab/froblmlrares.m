function f = froblmlrares(T, U, S, varargin)
%FROBLMLRARES Frobenius norm of residual of a LMLRA.
%   F = FROBLMLRARES(T, U) computes the Frobenius norm of the residual tensor E
%   which is given by lmlragen(U)-T.
%
%   FROBLMLRARES(T,U,options) or FROBLMLRARES(T,U,'key',value) can be used to set
%   the following options:
%
%   - ExpandLimit           Expands structured tensors to full tensors if the
%                           total number of entries prod(getsize(T)) is
%                           smaller than the value of this limit. Default:
%                           1e6.
%   - Type                  The type of the tensor. If not given, the type is
%                           determined automatically. 
%
%   All other options are passed on to LMLRARES if needed.
%
%   See also lmlrares, frob, getstructure
    
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
        f = frob(lmlrares(T, U, S, passopt));
    elseif prod(getsize(T)) < options.ExpandLimit
        if ~isfield(passopt, 'Format'), passopt.Format = false; end
        f = frob(lmlrares(ful(T),U,S,passopt));
    else 
        f = frob(T,'squared') - 2*real(inprod(T,{U,S})) + frob({U,S},'squared');
        f = sqrt(abs(f));
    end
end
