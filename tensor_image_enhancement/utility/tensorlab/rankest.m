function [R,L_lb,L_cpd] = rankest(T,varargin)
%RANKEST Estimate rank.
%   rankest(T) plots an L-curve of number of rank-one terms in a canonical
%   polyadic decomposition. The x-axis corresponds to the number of
%   rank-one terms, and the y-axis corresponds to the relative error of the
%   CPD in that many rank-one terms. Additionally, the corner R of the
%   resulting L-curve is estimated. The rank of a tensor can be different
%   when computed over the real field or the complex field. (See options
%   below.)
%
%   R = rankest(T) does not plot anything and instead returns the number of
%   rank-one terms R corresponding to the corner of the L-curve.
%
%   [R,L_lb,L_cpd] = rankest(T) also returns the L-curve ranks L_lb(:,1)
%   and L_cpd(:,1) and the corresponding lower bound on the truncation
%   error L_lb(:,2) and relative error of the CPD approximation L_cpd(:,2),
%   respectively.
%
%   rankest(T,options) may be used to set the following options:
%
%      options.MaxR =           - Maximum number of rank-one terms to try.
%      numel(T)/max(size_tens)
%      options.MinR = 1         - Minimum number of rank-one terms to try.
%                                 Only used if explicitly given, otherwise
%                                 MinR is chosen based on MaxRelErr, see
%                                 below. 
%      options.MaxRelErr = 1e-1 - Determines a lower threshold for the
%                                 number of rank-one terms R to try, based
%                                 on a lower bound of the relative error.
%      options.MinRelErr = 1e-2 - Determines an upper threshold for the
%                                 number of rank-one terms R to try.
%      options.Solver = @cpd    - The solver used to compute the CPD for
%                                 each number of rank-one terms R. Called
%                                 as options.Solver(T,R, ...
%                                 options.SolverOptions), where
%                                 options.SolverOptions is an options
%                                 structure passed to the solver.
%      options.XMultiplier = 1  - The importance of the number of rank-one
%                                 terms in determining the L-curve corner,
%                                 relative to the importance of the
%                                 relative error of the approximation.
%      options.Complex =        - Compute the rank over the real or the
%         [{'auto'}, true]        complex field. If 'auto', the rank is
%                                 computed over the real field if the tensor
%                                 is real-valued, and over the complex field
%                                 if the tensor is complex-valued. If equal
%                                 to true, the rank is computed over the
%                                 complex field.
%
%   See also mlrankest.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] J.L. Castellanos, S. Gomez, V. Guerra, "The triangle method for
%       finding the corner of the L-curve," Applied Numerical Mathematics,
%       Vol. 43, No. 4, 2002, pp. 359-373.
%
% Version History:
% - 2016/03/14   NV      Restored MaxRelErr option
% - 2015/10/14   NV      Added explicit option for complex tensors

type = getstructure(T);
if any(strcmpi(type, {'full', 'incomplete', 'sparse'}))
    T = fmt(T,true);
    type = getstructure(T);
end 
size_tens = getsize(T);

% Parse options
p = inputParser();
p.addOptional('MaxR', prod(size_tens)/max(size_tens));
p.addOptional('MinR', 1);
p.addOptional('MaxRelErr', 1e-1);
p.addOptional('MinRelErr', 1e-2);
p.addOptional('Complex', 'auto');
p.addOptional('Solver', @cpd);
p.addOptional('SolverOptions', {});
p.addOptional('XMultiplier', 1);
p.parse(varargin{:});
options = p.Results;

% Update some options
if isstruct(options.SolverOptions)
    fn = fieldnames(options.SolverOptions);
    fields = struct2cell(options.SolverOptions);
    tmp = [fn(:).'; fields(:).'];
    options.SolverOptions = tmp(:).';
end
fn = options.SolverOptions(1:2:end);
if ~any(strcmpi(fn,'ExploitStructure'))
    options.SolverOptions{end+1} = 'ExploitStructure';
    options.SolverOptions{end+1} = false;
end
if ~any(strcmpi(fn,'Display'))
    options.SolverOptions{end+1} = 'Display';
    options.SolverOptions{end+1} = false;
end
if ~any(strcmpi(fn,'Complex'))
    options.SolverOptions{end+1} = 'Complex';
    options.SolverOptions{end+1} = options.Complex;
end
if ~any(strcmpi(fn,'TolFun'))
    if options.MinRelErr < 1e-6, 
        TolFun = max(options.MinRelErr^2, eps^2);
        if options.MinRelErr < 1e-14
            TolFun = eps^2;
        end
        options.SolverOptions{end+1} = 'TolFun';
        options.SolverOptions{end+1} = TolFun;
    end
end
if ~any(strcmpi(fn,'TolX'))
    if options.MinRelErr < 1e-6, 
        TolX = max(options.MinRelErr, eps/10);
        if options.MinRelErr < 1e-13
            TolX = 0;
        end
        options.SolverOptions{end+1} = 'TolX';
        options.SolverOptions{end+1} = TolX;
    end
end

% Compute lower bound on truncation error.
frobT = frob(T);
if any(strcmpi(type, {'full', 'sparse'}))
    if isstruct(T)
        [~,~,sv] = mlsvds(T);
    else 
        [~,~,sv] = mlsvd(T);
    end
    Rmax = max(cellfun(@length,sv))-1;
    lb = max(cell2mat(cellfun( ...
             @(s)[sqrt(flipud(cumsum(flipud(s(:).^2)))).'/frobT ...
             zeros(1,Rmax+1-length(s))],sv(:),'UniformOutput',false)));
    lb = lb(2:end);
else
    Rmax = 0;
    lb = [];
end

% Determine the range of rank-one terms to test.
if ~isempty(lb) && any(strcmpi(p.UsingDefaults, 'MinR'))
    %R = find(lb(options.MinR:end)<=options.MinRelErr,1,'first')+ ...
    R = find(lb(options.MinR:end)<=options.MaxRelErr,1,'first') + ...
        options.MinR-1;
    if isempty(R), R = max(options.MinR,Rmax+1); end
    R = R:options.MaxR;
else
    R = 1:options.MaxR;
end

% Compute the relative error for each number of rank-one terms R(r).
relerr = zeros(1,length(R));
warningstate = warning('query', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix')
for r = 1:length(R)
    
    % Compute the CPD in R rank-one terms.    
    [U,output] = options.Solver(T,R(r),options.SolverOptions{:});
    if isfield(output,'Refinement') && ...
            isfield(output.Refinement, 'relerr')
        relerr(r) = output.Refinement.relerr;
    elseif isfield(output,'Refinement') && ...
       isfield(output.Refinement,'fval')
        relerr(r) = sqrt(2*output.Refinement.fval(end))/frobT;
    elseif isfield(output,'Algorithm') && ...
            isfield(output.Algorithm, 'relerr')
        relerr(r) = output.Algorithm.relerr;
    elseif isfield(output,'Algorithm') && ...
           isfield(output.Algorithm,'fval')
        relerr(r) = sqrt(2*output.Algorithm.fval(end))/frobT;
    elseif isfield(output,'fval')
        relerr(r) = sqrt(2*output.fval(end))/frobT;
    else
        relerr(r) = frobcpdres(T,U)/frobT;
    end
    
    % Stop is options.MinRelErr has been reached.
    if relerr(r) <= options.MinRelErr
        R = R(1:r);
        relerr = relerr(1:r);
        break;
    end
    
    % Update the plot.
    if r > 1 && r < length(R) && nargout == 0, 
        if r == 2, clf; end
        Lcurve(true); 
    end

end
warning(warningstate.state, 'MATLAB:nearlySingularMatrix')

% Compute L-curve corner using the triangle method [1].
logx = options.XMultiplier*[R(:);R(end)+1];
logy = [log10(relerr(:));log10(relerr(end))];
if R(1) ~= 1
    logx = [options.XMultiplier;logx];
    if ~isempty(lb), logy = [log10(lb(1));logy];
    else logy = [0;logy]; end
end
ab = cat(3,bsxfun(@minus,logx,logx.'),bsxfun(@minus,logy,logy.'));
ac = cat(3,logx(end)-logx.',logy(end)-logy.');
area = bsxfun(@times,ab(:,:,1),ac(:,:,2)) - ...
       bsxfun(@times,ab(:,:,2),ac(:,:,1));
cosa = bsxfun(@times,ab(:,:,1),ac(:,:,1)) + ...
       bsxfun(@times,ab(:,:,2),ac(:,:,2));
cosa = bsxfun(@rdivide,cosa./sqrt(ab(:,:,1).^2+ab(:,:,2).^2), ...
                       sqrt(ac(:,:,1).^2+ac(:,:,2).^2));
cosa(area >= 0 | tril(true(size(cosa))) | cosa <= cos(7*pi/8)) = -1;
[a,opt] = max(max(cosa));
if isempty(a) || a == -1, opt = 1;
elseif R(1) ~= 1, opt = opt-1; end

% Display output.
L_lb = [(1:Rmax).' lb(:)];
L_cpd = [R(:) relerr(:)];
if nargout == 0, Lcurve(false); end
R = R(opt);

function Lcurve(update)
    if R(r) <= 1, return; end
    style = {'Marker','+','MarkerSize',2.5};
    if ~update
        semilogy(R(opt),relerr(opt),'rs','LineStyle','none'); hold on;
    end
    semilogy(R(1:length(relerr)),relerr,'r',style{:}); hold on;
    semilogy(1:length(lb),lb,style{:});
    semilogy([1 R(r)],options.MaxRelErr*[1 1],'k:');
    semilogy([1 R(r)],options.MinRelErr*[1 1],'k:');
    text(1,options.MaxRelErr,'MaxRelErr','VerticalAlignment','Top');
    text(1,options.MinRelErr,'MinRelErr','VerticalAlignment','Top');
    hold off;
    if isempty(lb), lb1 = options.MaxRelErr*10^0.5; else lb1 = lb(1); end
    ylim([min(options.MinRelErr/(10^0.5),min(relerr(1:r))) ...
          max(options.MaxRelErr*10^0.5,lb1)]);
    %ylim([min(options.MinRelErr/(10^0.5),min(relerr(1:r))) 1]);
    xlim([1 R(r)]);
    ylabel('frob(cpdres(T,U))/frob(T)'); xlabel('R');
    xt = get(gca,'XTick'); set(gca,'XTick',xt(mod(xt,1) == 0));
    if update
        if ~isempty(lb)
            legend(['CPD error (trying R = ' int2str(R(r+1)) '...)'], ...
                'Lower bound on error','Location','NE');
        else
            legend(['CPD error (trying R = ' int2str(R(r+1)) '...)'], ...
                'Location','NE');
        end
    else
        if ~isempty(lb)
            legend(['L-curve corner at R = ' int2str(R(opt))], ...
                'CPD error','Lower bound on error','Location','NE');
        else
            legend(['L-curve corner at R = ' int2str(R(opt))], ...
                'CPD error','Location','NE');
        end
    end
    set(datacursormode(gcf),'UpdateFcn',@datacursor);
    drawnow;
end

function txt = datacursor(~,event_obj)
    pos = get(event_obj,'Position');
    txt = {['R: ' int2str(pos(1))], ...
           ['relative error: ' num2str(pos(2))]};
end

end
