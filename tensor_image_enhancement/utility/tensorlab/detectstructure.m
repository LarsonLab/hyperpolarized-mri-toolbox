function [structures,representations] = detectstructure(T,varargin)
%DETECTSTRUCTURE Detects structure in a tensor
%   structures = DETECTSTRUCTURE(T) detects (block-)Hankel, symmetry and conjugate
%   symmetry structure in the tensor T. 
%
%   [structures,representations] = DETECTSTRUCTURE(T) also returns the
%   efficient representations of the structured versions, if structure has
%   been detected.
%
%   DETECTSTRUCTURE(T,'key',value,...) or DETECTSTRUCTURE(T,options) can be
%   used to pass the following options:
%
%   - Structures:    A string, or cell of strings, indicating the
%                    structure(s) the method searches for. Currently
%                    supported:
%                       - 'hankel' for Hankel structure
%                       - 'symmetric' for symmetry
%                       - 'csymmetric' for conjugate symmetry
%
%   - RelErr:        The threshold for the relative error in each equality
%                    imposed by the structure. All relative errors should
%                    be lower than the threshold such that the structure is
%                    detected. Default: 1e-13.

%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version History:
%   - 2016/01/10   OD      Initial version

supportedstructures = {'hankel','symmetric','csymmetric'};

p = inputParser;
p.addOptional('Structures', supportedstructures);
p.addOptional('RelErr',1e-13);
p.parse(varargin{:});
options = p.Results;

type = getstructure(T);
if ~strcmp(type,'full')
    error('detectstructure:T','Non-full tensors are not supported yet!');
end

% Convert options.Structures to cell
wascell = true;
if ~iscell(options.Structures), wascell = false; options.Structures = {options.Structures}; end
options.Structures = cellfun(@lower,options.Structures,'UniformOutput',false);

structures = {};
representations = {};

%% detect Hankel structure
Tidxn0 = T~=0;
msizeT = cell(1,ndims(T));
for t = 1:ndims(T)
    if size(T,t)<3, msizeT{t} = ':';
    else msizeT{t} = 1:3;
    end
end
Tidxn0small = Tidxn0(msizeT{:});

% Measure on the small part of the data
measuresmall = @(H,Hh) all(abs(H(Tidxn0small)-Hh(Tidxn0small))<=options.RelErr*abs(H(Tidxn0small))) && ...
    all(abs(H(~Tidxn0small)-Hh(~Tidxn0small))<=options.RelErr);
% Measure on the large part of the data
measurelarge = @(H,Hh) all(abs(H(Tidxn0)-Hh(Tidxn0))<=options.RelErr*abs(H(Tidxn0))) && ...
    all(abs(H(~Tidxn0)-Hh(~Tidxn0))<=options.RelErr);

if any(strcmp(options.Structures,'hankel'))
    % Check for Hankel tensor structure
    [ishankelt,repr] = ishankeltensor(T);
    if ishankelt
        structures{end+1} = 'hankel';
        representations{end+1} = repr;
    else
        % Check for block-Hankel structure
        [ishankelb,repr] = isanyblockhankeltensor(T);
        if ishankelb
            structures{end+1} = 'hankel';
            representations{end+1} = repr;
        end
    end
end

%% detect symmetry
if any(strcmp(options.Structures,'symmetric'))
    if ~isempty(structures) && strcmp(structures{end},'hankel')
        structures{end+1} = 'symmetric';
        representations{end+1} = [];
    else
        if ismatrix(T)
            if size(T,1)==size(T,2) && frob(T.'/2-T/2)/frob(T)<options.RelErr
                structures{end+1} = 'symmetric';
                representations{end+1} = [];
            end
        end
    end
end

%% detect conjugate symmetry
if any(strcmp(options.Structures,'csymmetric'))
    if ismatrix(T)
        if size(T,1)==size(T,2) && frob(T'/2-T/2)/frob(T)<options.RelErr
            structures{end+1} = 'csymmetric';
            representations{end+1} = [];
        end
    end
end

if ~wascell && ~isempty(structures)
    structures = structures{1};
    representations = representations{1};
end

%% Auxiliary functions
    function [flag,representation] = isanyblockhankeltensor(H)
        % Detect if H is a block-Hankel tensor
        
        flag = false;
        representation = [];
        N = ndims(H);
        for i = N-1:-1:2
            % Iterate from high to low Hankel order
            modes = nchoosek(1:N,i);
            for j = 1:size(modes,1)
                [ismostblockhankel,representation] = ismostblockhankeltensor(H,modes(j,:));
                if ismostblockhankel
                    % block-Hankel structure has been detected
                    flag = true;
                    repermorder = [];
                    repermorder(modes(j,:)) = 1:i;
                    othermodes = 1:N; othermodes(modes(j,:)) = [];
                    repermorder(othermodes) = i+1:N;
                    representation.repermorder = repermorder;
                    return;
                end
            end
        end
    end

    function [flag,representation] = ismostblockhankeltensor(H,dims)
        % Detect if H has Hankel structure in the given dimensions dims
        
        flag = false;
        representation = [];
        N = ndims(H);
        
        % Check small
        msize = cell(1,ndims(H));
        for i = 1:ndims(H)
            if size(H,i)<3, msize{i} = ':';
            else msize{i} = 1:3;
            end
        end
        tmp = ful(H,msize{:});
        
        X = dehankelize(tmp,'method','fibers','dims',dims);
        sz = size(tmp);
        ind = cumsum(sz(dims))-(0:(numel(dims)-1));
        Hh = hankelize(X,'dim',dims(1),'ind',ind(1:end-1));
        otherdims = 1:N; otherdims(dims) = [];
        if dims(1)>numel(otherdims)
            permorder = [1:(dims(1)-1) dims];
        else
            permorder = [1:(dims(1)-1) dims otherdims(dims(1):end)];
        end
        
        tmpperm = permute(tmp,permorder);
        if ~measuresmall(tmpperm,Hh)
            return;
        end
        
        % Check large
        X = dehankelize(H,'method','fibers','dims',dims);
        sz = size(H);
        ind = cumsum(sz(dims))-(0:(numel(dims)-1));
        [Hh,Hhs] = hankelize(X,'dim',dims(1),'ind',ind(1:end-1));
        H = permute(H,permorder);
        flag = measurelarge(H,Hh);
        representation = Hhs;
        
    end

    function [flag,representation] = ishankeltensor(H)
        % Detect if H is a Hankel tensor in all of the modes
        
        flag = false;
        representation = [];
        
        % Check small
        msize = cell(1,ndims(H));
        for i = 1:ndims(H)
            if size(H,i)<3, msize{i} = ':';
            else msize{i} = 1:3;
            end
        end
        tmp = ful(H,msize{:});
        X = dehankelize(tmp,'method','fibers','order',ndims(H));
        sz = size(tmp);
        ind = cumsum(sz(1:end-1))-(0:(ndims(tmp)-2));
        Hh = hankelize(X,'ind',ind);
        
        if ~measuresmall(tmp,Hh)
            return;
        end
        
        % Check large
        X = dehankelize(H,'method','fibers','order',ndims(H));
        sz = size(H);
        ind = cumsum(sz(1:end-1))-(0:(ndims(H)-2));
        [Hh,Hhs] = hankelize(X,'ind',ind);
        flag = measurelarge(H,Hh);
        representation = Hhs;
    end

end