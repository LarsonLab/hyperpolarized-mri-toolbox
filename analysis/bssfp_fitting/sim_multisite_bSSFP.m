
function [Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, varargin)
% Inputs:
% flips - struct with vector members P, L, B, A etc representing flip angles in
%     degrees for each excitation per time point, e.g. flips.P should be 
%     size [Nt num_ex]
%     NOTE: for 2D acquisition num_ex=1, i.e. flips is of size [Nt 1]
% TR - strcut with vector members P, L, B, A etc representing time between data
%     excitations for each excitation in seconds e.g. TR.P should be size [1 num_ex]
%     NOTE: for 2D acquisition num_ex=1, i.e. a scalar and TR should be
%     equivalent to total time spent acquiring full metabolite image
%     Also, sum of TRs across all metabolites ~= TempRes
% TempRes - scalar representing total time between time points in seconds
% R1 - vector with size [1 Nmets] with the reciprocal of T1 in seconds 
%     for each metabolite
% k - vector for forward apparent rate constant- size [nmets-1 1] for ex.
%     for two site should be value for [kPL] in 1/s
% Nt - # of time points, scalar
% Mz0 - starting longitudinal magnetization for each metabolite,
%       vector [1 Nmets]
%
% Optional Inputs:
% R2 - struct with sclar members P, L, B, A representing the reciprocal of T2 
%   in s for pyruvate/lactate, NOTE: should be included for bSSFP sequence
% spoilers - vector with with values equal to 1 or 0 
%   specifying whether a spoiler was used for each metabolite in order P, L, A, B, 
%   0=no spoiler (bSSFP), 1=spoiler (GRE) (default: 1) 
% scales - vector of values to use to scale the signal after modelling of
%   size [1 Nmets] (default: [1 1 1 1])
% cat_flips - vector for catalyzation flip angle sequence for bSSFP 
%    acquisition, default is [half of 1st bSSFP flip angle]
% cat_TR - vector for catalyzation TR sequence for bSSFP acquisition,
%    default is [half of bSSFP TR] 
% input - input function [1 Nt], if input isnt provided will assume no
%   input
% plot_flag - 1 to display plots, 0 to not, default:0
%
% Outputs:
% Mxy - a complex array of [Nmets Nt] signal/ transverse mag
% Mz - an array of [Nmets Nt] longitudinal mag before RF
%
% Author: Sule Sahin, 2024 Copyright

    % TO DO: add checks in place to catch errors
    p = inputParser;
    p.addParameter( 'R2', [], @isstruct );
    p.addParameter( 'spoilers', [], @isvector );
    p.addParameter( 'scales', [], @isvector );
    p.addParameter( 'cat_flips', [] );
    p.addParameter( 'cat_TR', [] );
    p.addParameter( 'input', [] );
    p.addParameter( 'plot_flag', false, @(x) islogical(x) || isnumeric(x) );
    p.parse( varargin{:} );
    R2 = p.Results.R2;
    spoilers = p.Results.spoilers;
    scales = p.Results.scales;
    cat_flips = p.Results.cat_flips;
    cat_TR = p.Results.cat_TR;
    input = p.Results.input;
    plot_flag = p.Results.plot_flag;
    
    nmets = numel(fieldnames(flips));
    
    if ~isstruct(R2) %if no R2 given set it equal to 0
        R2.P = 0; R2.L = 0; R2.B = 0; R2.A=0;
    end
    if ~sum(strcmp(fieldnames(R2), 'P'))
        R2.P = 0;
    end
    if ~sum(strcmp(fieldnames(R2), 'L'))
        R2.L = 0;
    end
    if ~sum(strcmp(fieldnames(R2), 'B'))
        R2.B = 0;
    end
    if ~sum(strcmp(fieldnames(R2), 'A'))
        R2.A = 0;
    end
        
    if isempty(spoilers) %if no spoilers given, assume spoiled sequence for all mets
        spoilers = [1 1 1 1];
    end

    if isempty(scales) %if no scales given, assume no scaling in b/w mets
        scales = ones([1 nmets]);
    end
    
    % catalyzation pulse set up
    cat_flips_tot = cell(1,nmets);
    cat_TR_tot = cell(1,nmets);
    if ~all(spoilers) %if any bSSFP sequence
        met_list={'P', 'L', 'B', 'A'};
        ssfp_idx = find(spoilers==0); %find # PE for bssfp using TR
        flips_ssfp = flips.(met_list{ssfp_idx(1)});
        TR_ssfp = TR.(met_list{ssfp_idx(1)});
        for ii=1:length(ssfp_idx)
            if or(isempty(cat_flips), isempty(cat_TR)) %if either empty set to defaults
                cat_flips_tot{ssfp_idx(ii)} = [flips_ssfp(1)/2];
                cat_TR_tot{ssfp_idx(ii)} = [TR_ssfp(1)/2];
            else
                cat_flips_tot{ssfp_idx(ii)} = cat_flips;
                cat_TR_tot{ssfp_idx(ii)} = cat_TR;
            end
        end
    end

    
    if size(flips.P, 1) == 1
        flips.P = repmat(flips.P, [Nt, 1]);
    end
    if size(flips.L, 1) == 1
        flips.L = repmat(flips.L, [Nt, 1]);
    end
    if isfield(flips, 'B') && size(flips.B, 1) == 1
        flips.B = repmat(flips.B, [Nt, 1]);
    end
    if isfield(flips, 'A') && size(flips.A, 1) == 1
        flips.A = repmat(flips.A, [Nt, 1]);
    end
    
    if isempty(input)  % TO DO : cleanup
        input = zeros([1 Nt]);
    end
    input_fxn = zeros([nmets*3 Nt]);    
    if size(input) == [nmets Nt]
        input_fxn(3,:) = input(1,:);
        input_fxn(6,:) = input(2,:);
        if nmets>2
            input_fxn(9,:) = input(3,:);
        end
        if nmets>3
            input_fxn(12,:) = input(4,:);
        end
    elseif size(input) == [nmets*3 Nt]
        input_fxn = input;
    else       
        input_fxn(3, :) = input;
    end

    %check TR and TempRes defined correctly
    TotalTR = 0;
    TRfields = fieldnames(TR);
    for f=1:numel(TRfields)
        TotalTR = TotalTR + sum(TR.(TRfields{f}));
    end
    TotalTR = TotalTR + sum(cat(2,cat_TR_tot{:}))*2;
    if abs(TempRes - TotalTR) > 0.3
        warning("Sum of TRs is not equal or close to Temporal Resolution. The difference between the sum and Temp Res is more than 0.3s. Please check definition of TR and TempRes.")
    end
    
%     % determine FA scaling
%     FA_scale = ones([1 nmets]);
%     met_list={'P', 'L', 'B', 'A'};
%     if any(spoilers) %if any gre sequence (if not dont need to scale)
%         if all(spoilers)
%             Nb = size(flips.P,2); %if only GRE use # pe for pyruvate
%         else
%             ssfp_idx = find(spoilers==0); %find # PE for bssfp using TR
%             TR_ssfp = TR.(met_list{ssfp_idx(1)});
%             TR_nonramp = TR_ssfp(round(length(TR_ssfp)/2)); %dont include catalyzation pulses
%             Nb = size(TR_ssfp(TR_ssfp==TR_nonramp),2); 
%         end
%         for met= 1:nmets
%             if spoilers(met) == 1
%                 flips_met = flips.(met_list{met});
%                 Npe = size(flips_met,2);
%                 FA = flips_met(1,1);
%                 FA_scale(met) = sind(acosd(nthroot(cosd(FA)^Npe,Nb))) / sind(FA);
%             end
%         end
%     end

    % determine bssfp vs GRE scaling
%     Tread_scale = ones([1 nmets]);
%     if any(spoilers) %if any gre sequence (if not dont need to scale)
%         gre_idx = find(spoilers); %scale by Tread of 1st met w/ gre sequence
%         for met= 1:nmets
%             if spoilers(met) == 0
%                 %Tread_scale(met) = Tread(met) ./ Tread(gre_idx(1)); %original scale
%                 Tread_scale(met) = Tread(gre_idx(1)) ./ Tread(met); %inverse/theoretical scale
%             end
%         end
%     end
    
    %define R for given relaxation rates and initial mag vector
    switch nmets
        case 2
            M0 = [0; 0; Mz0(1); 0; 0; Mz0(2)]; 
            R = [-R2.P  0      0       0    0    0;
                   0  -R2.P    0       0    0    0;
                   0    0  -k(1)-R1(1)  0    0    0;
                   0    0      0     -R2.L  0    0;
                   0    0      0       0  -R2.L  0;
                   0    0      k(1)    0    0  -R1(2)];
            met_list = {'pyruvate', 'lactate'};
        case 3
            M0 = [0; 0; Mz0(1); 0; 0; Mz0(2); 0; 0; Mz0(3)]; 
            R = [-R2.P  0        0          0    0    0    0    0    0;
                   0  -R2.P      0          0    0    0    0    0    0;
                   0    0  -k(1)-k(2)-R1(1)  0    0    0    0    0    0;
                   0    0        0        -R2.L  0    0    0    0    0;
                   0    0        0          0  -R2.L  0    0    0    0;
                   0    0        k(1)       0    0  -R1(2)  0    0    0;
                   0    0        0          0    0    0  -R2.B 0    0;
                   0    0        0          0    0    0    0  -R2.B  0;
                   0    0        k(2)       0    0    0    0    0  -R1(3)];
             met_list = {'pyruvate', 'lactate', 'bicarb'};
         case 4
            M0 = [0; 0; Mz0(1); 0; 0; Mz0(2); 0; 0; Mz0(3); 0; 0; Mz0(4)]; 
            R = [-R2.P  0            0           0    0    0    0    0    0    0    0    0;
                   0  -R2.P          0           0    0    0    0    0    0    0    0    0;
                   0    0  -k(1)-k(2)-k(3)-R1(1)  0    0    0    0    0    0    0    0    0;
                   0    0            0         -R2.L  0    0    0    0    0    0    0    0;
                   0    0            0           0  -R2.L  0    0    0    0    0    0    0;
                   0    0            k(1)        0    0  -R1(2)  0    0    0    0    0    0;
                   0    0            0           0    0    0  -R2.B  0    0    0    0    0;
                   0    0            0           0    0    0    0  -R2.B  0    0    0    0;
                   0    0            k(2)        0    0    0    0    0  -R1(3)  0    0    0;
                   0    0            0           0    0    0    0    0    0  -R2.A  0    0;
                   0    0            0           0    0    0    0    0    0    0  -R2.A  0;
                   0    0            k(3)        0    0    0    0    0    0    0    0  -R1(4)];
            met_list = {'pyruvate', 'lactate', 'bicarb', 'alanine'};
    end

    %propogate magnetization vector
    M_plus = M0;
    Mt = zeros([size(M0,1) Nt]);
    for t=1:Nt
        M_plus = M_plus + input_fxn(:,t);   

        %pyruvate
        if numel(flips.P(t,:)) > 1
            [Mt(1:3, t), M_next] = propogate_train_3D(flips.P(t,:), R, TR.P, spoilers(1), 'P', nmets, scales(1), cat_flips_tot{1}, cat_TR_tot{1}, M_plus);
        else
            [Mt(1:3, t), M_next] = propogate_train_2DGRE(flips.P(t,:), R, TR.P, spoilers(1), 'P', nmets, scales(1), M_plus);
        end
        
        
        %lactate
        if numel(flips.L(t,:)) > 1
            [Mt(4:6,t), M_next] = propogate_train_3D(flips.L(t,:), R, TR.L, spoilers(2), 'L', nmets, scales(2), cat_flips_tot{2}, cat_TR_tot{2}, M_next);
        else
            [Mt(4:6,t), M_next] = propogate_train_2DGRE(flips.L(t,:), R, TR.L, spoilers(2), 'L', nmets, scales(2), M_next);
        end
        
        %bicarb
        if nmets > 2
            if numel(flips.B(t,:)) > 1
                [Mt(7:9,t), M_next] = propogate_train_3D(flips.B(t,:), R, TR.B, spoilers(3), 'B', nmets, scales(3), cat_flips_tot{3}, cat_TR_tot{3}, M_next);
            else
                [Mt(7:9,t), M_next] = propogate_train_2DGRE(flips.B(t,:), R, TR.B, spoilers(3), 'B', nmets, scales(3), M_next);
            end
        end
        
        %alanine
        if nmets > 3
            if numel(flips.A(t,:)) > 1
                [Mt(10:12,t), M_next] = propogate_train_3D(flips.A(t,:), R, TR.A, spoilers(4), 'A', nmets, scales(4), cat_flips_tot{4}, cat_TR_tot{4}, M_next);
            else
                [Mt(10:12,t), M_next] = propogate_train_2DGRE(flips.A(t,:), R, TR.A, spoilers(4), 'A', nmets, scales(4), M_next);
            end
        end
        
        M_plus = M_next;

    end

    Mxy = abs(Mt(1:3:end, :)) + 1i.*abs(Mt(2:3:end, :));
    Mz = Mt(3:3:end,:);

    % plot fits vs measured data
    if plot_flag
        t=0:TempRes:(Nt-1)*TempRes;
        figure('Position',[300,500,500,500]);
        subplot(211)
        plot(t, abs(Mxy), 'LineWidth', 2)
        title('|Mxy|')
        xlabel('time (s)')
        ylabel('|Mxy|')
        legend(met_list)
        
        subplot(212)
        plot(t, Mz, 'LineWidth', 2)
        title('Mz')
        xlabel('time (s)')
        ylabel('Mz')
        legend(met_list)
    end
end

function [M_t, M_next] = propogate_train_2DGRE(flips, R, TR, spoilers, met, nmets, norm, M_next)
    
    spoil = diag(repmat([0 0 1], [1 nmets]))^spoilers;
    M_ti = rot_matrix_met(flips, met, nmets) * spoil * M_next;
    switch met
        case 'P'
            M_t = [norm .* M_ti(1:2); M_next(3)];
        case 'L'
            M_t = [norm .* M_ti(4:5); M_next(6)];
        case 'B'
            M_t = [norm .* M_ti(7:8); M_next(9)];
        case 'A'
            M_t = [norm .* M_ti(10:11); M_next(12)];
    end
    M_next = spoil * expm(R*TR) * M_ti;
end


function [M_t, M_next] = propogate_train_3D(flips_nocat, R, TR_nocat, spoilers, met, nmets, norm, cat_flips, cat_TR, M_next)
    
    spoil = diag(repmat([0 0 1], [1 nmets]))^spoilers;
    cat_len = size(cat_flips,2);
    num_ex = size(flips_nocat,2) + 2*cat_len;
    start_idx = cat_len + 1;
    end_idx = num_ex - cat_len;
    
    if rem(size(flips_nocat,2),2) %odd number of flips
        flips = cat(2, cat_flips, flips_nocat, fliplr(cat_flips));
    else %even number
        flips = cat(2, cat_flips, flips_nocat, fliplr(cat_flips) .* -1);
    end
    TR = cat(2, cat_TR, TR_nocat, fliplr(cat_TR));
    
    M_ti = zeros([3 length(flips_nocat)]);
    for i=1:num_ex
        rot = rot_matrix_met(flips(i), met, nmets);
        if i==1
            Mz = M_next;
            M_next = rot * spoil * M_next;
        else
            M_next = rot * expm(R*TR(i-1)) * spoil * M_next;
        end
        if and(i >= start_idx, i <= end_idx)
            switch met
                case 'P'
                    M_ti(:,i-cat_len) = M_next(1:3);
                case 'L'
                    M_ti(:,i-cat_len) = M_next(4:6);
                case 'B'
                    M_ti(:,i-cat_len) = M_next(7:9);
                case 'A'
                    M_ti(:,i-cat_len) = M_next(10:12);
            end
        end
    end
    M_next = expm(R*TR(end)) * spoil * M_next; %last bit of relaxation after last flip
    switch met
    case 'P'
        Mz = Mz(3);
    case 'L'
        Mz = Mz(6);
    case 'B'
        Mz = Mz(9);
    case 'A'
        Mz = Mz(12);
    end

    M_t = [norm .* mean(abs(M_ti(1:2,:)),2); Mz]; 
end


function rot = rot_matrix_met(FA, met, nmets)
    switch met
        case 'P'
            rot = [[rot_matrix(FA) zeros(3) zeros(3) zeros(3)]; 
                   [zeros(3)       eye(3)   zeros(3) zeros(3)];
                   [zeros(3)       zeros(3) eye(3)   zeros(3)];
                   [zeros(3)       zeros(3) zeros(3) eye(3)]];
        case 'L'
            rot = [[eye(3)   zeros(3)       zeros(3) zeros(3)]; 
                   [zeros(3) rot_matrix(FA) zeros(3) zeros(3)];
                   [zeros(3) zeros(3)       eye(3)   zeros(3)];
                   [zeros(3) zeros(3)       zeros(3) eye(3)]];
        case 'B'
            rot = [[eye(3)   zeros(3) zeros(3)       zeros(3)]; 
                   [zeros(3) eye(3)   zeros(3)       zeros(3)];
                   [zeros(3) zeros(3) rot_matrix(FA) zeros(3)];
                   [zeros(3) zeros(3) zeros(3)       eye(3)]];
        case 'A'
            rot = [[eye(3)   zeros(3) zeros(3) zeros(3)]; 
                   [zeros(3) eye(3)   zeros(3) zeros(3)];
                   [zeros(3) zeros(3) eye(3)   zeros(3)];
                   [zeros(3) zeros(3) zeros(3) rot_matrix(FA)]];
    end
    if nmets==2
        rot = rot(1:6, 1:6);
    elseif nmets==3
        rot = rot(1:9, 1:9);
    end
end

function rot = rot_matrix(FA)

    rot = [ cosd(FA) 0 -sind(FA);
        0 1 0;
        sind(FA) 0 cosd(FA)];

end