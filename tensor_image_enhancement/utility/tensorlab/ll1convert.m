function [Uout,L] = ll1convert(Uin,varargin)
%LL1CONVERT Convert LL1 decomposition between CPD and BTD format.
%   UCPD = LL1CONVERT(UBTD) converts a decomposition in LL1 terms given in the
%   BTD format to the CPD format. The BTD format UBTD is a cell with R =
%   length(L) terms, each having two factor matrices with L(r) columns, a column
%   vector, and a L(r)xL(r) matrix. The CPD format UCPD is a cell with 3
%   factor matrices, the first two having sum(L) columns, the third one
%   having R=length(L) columns.
%
%   [UCPD,L] = LL1CONVERT(UBTD) also returns the computed L.
%
%   UBTD = LL1CONVERT(UCPD, L) converts a decomposition in LL1 terms given in
%   the CPD format to the BTD format.   
%
%   LL1CONVERT(U,'OutputFormat',format) or LL1CONVERT(U,L,'OutputFormat',
%   format) can be used to explicitly set the output format to 'btd' or
%   'cpd'.
    
%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
    
    if nargin >= 2 && isnumeric(varargin{1})
        L = varargin{1}(:).'; 
        varargin = varargin(2:end);
    end
    p = inputParser;
    p.addOptional('OutputFormat', 'auto');
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    Uin = Uin(:).';
    
    if all(cellfun(@isnumeric,Uin)) % CPD format
        if nargin == 1
            error('ll1convert:L', ...
                  'L should be given if Ucpd is in the CPD format.');
        end 
        if length(Uin) ~= 3
            error('ll1convert:Ucpd', ...
                  ['Only third-order tensors can be handled, i.e, length(UCPD) ' ...
                   'should be 3.']);
        end
        if any(~cellfun(@ismatrix, Uin))
            error('ll1convert:Ucpd', ...
                  'Ucpd{n} should be matrices/vectors for all n.');
        end
        if ~all(cellfun('size', Uin, 2) == [sum(L), sum(L), length(L)])
            error('ll1convert:Ucpd', ['size(Ucpd{n},2) should be sum(L) for n=1,2 ' ...
                                'and length(L) for n=3.']);
        end
        
        inputFormat = 'cpd';
        if strcmpi(options.OutputFormat, 'auto')
            options.OutputFormat = 'btd';
        end
    elseif all(cellfun(@iscell, Uin)) % BTD format
        Uin = cellfun(@(u) u(:).', Uin, 'UniformOutput', false);
        if any(cellfun(@length, Uin) ~= 4)
            error('ll1convert:Ubtd', ['All terms Ubtd{r} should define third-order ' ...
                                'tensors, i.e., length(Ubtd{r}) should be 4 ' ...
                                'for all r. ']);
        end

        if exist('L', 'var')
            if length(Uin) ~= length(L)
                error('ll1convert:L', ['The detected number of terms length(Ubtd) ' ...
                                    'should match R=length(L).']);
            end
            if any(cellfun(@(u) size(u{1},2), Uin) ~= L)
                error('ll1convert:L', ['The detected L=cellfun(@(u) size(u{end},1), ' ...
                                'Ubtd) and the given L are not equal. (L ' ...
                                'is optional in the BTD format.)']);
            end 
        else
            L = cellfun(@(u) size(u{1},2), Uin);
        end
        if any(cellfun(@(u) any(~cellfun(@ismatrix, u)), Uin))
            error('ll1convert:Ubtd', ['All entries Ubtd{r}{n} should be matrices ' ...
                                'or vectors.']);
        end
        if any(cellfun(@(u) size(u{2}, 2), Uin) ~= L)
            error('ll1convert:Ubtd', ['All terms Ubtd{r}{2} should have L(r) ' ...
                                'columns, i.e., size(Ubtd{r}{2},2) should be ' ...
                                'size(Ubtd{r}{1},2).']);
        end
        if any(cellfun(@(u) size(u{3}, 2), Uin) ~= 1)
            error('ll1convert:Ubtd', ['All terms Ubtd{r}{3} should be column ' ...
                                'vectors']);
        end
        if any(cellfun(@(u) size(u{4}, 1), Uin) ~= L) || ...
                any(cellfun(@(u) size(u{4}, 2), Uin) ~= L)
            error('ll1convert:Ubtd', ['All terms U{r}{4} should be L(r)xL(r) ' ...
                                'matrices, in which L(r)=size(Ubtd{r}{1},2)']);
        end
        size_tens = cellfun('size', Uin{1}(1:3), 1);
        if any(cellfun(@(u) any(cellfun('size',u(1:3),1)~=size_tens), Uin))
            error('ll1convert:Ubtd', ['All terms Ubtd{r} should define tensors ' ...
                                'of the same size = cellfun(''size'',' ...
                                'Ubtd{1}(1:3),1).']);
        end
        
        inputFormat = 'btd';
        if strcmpi(options.OutputFormat, 'auto')
            options.OutputFormat = 'cpd';
        end
    else 
        error('ll1convert:U', 'Unknown input format');
    end
    
    if strcmpi(options.OutputFormat, 'btd')
        if strcmpi(inputFormat, 'cpd')
            % Do actual conversion
            Uin{1} = mat2cell(Uin{1}, size(Uin{1},1), L);
            Uin{2} = mat2cell(Uin{2}, size(Uin{2},1), L);
            Uin{3} = num2cell(Uin{3}, 1);
            Uout = cell(1,length(L));
            for r = 1:length(L)
                Uout{r} = {Uin{1}{r}, Uin{2}{r}, Uin{3}{r}, eye(L(r))};
            end
        else 
            Uout = Uin;
            for r = 1:length(Uin), 
                Uout{r}{1} = Uout{r}{1}*Uout{r}{4}; 
                Uout{r}{4} = eye(size(Uout{r}{4}));
            end
        end
    elseif strcmpi(options.OutputFormat, 'cpd')
        if strcmpi(inputFormat, 'btd')
            % Do actual conversion
            for r = 1:length(Uin), Uin{r}{1} = Uin{r}{1}*Uin{r}{4}; end
            Uout = cell(1,3);
            for n = 1:3
                Uout{n} = cellfun(@(u) u{n}, Uin, 'UniformOutput', false);
                Uout{n} = cat(2, Uout{n}{:});
            end
        else 
            Uout = Uin;
        end
    else 
        error('ll1convert:OutputFormat', 'Unknown output format %s.', ...
              options.OutputFormat);
    end
end
