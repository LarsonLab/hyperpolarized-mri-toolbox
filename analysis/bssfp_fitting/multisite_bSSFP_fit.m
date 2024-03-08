
function [finalfit, metrics, Mxy_fit, Mz_fit, input_all] = multisite_bSSFP_fit(S, params_fixed, params_est, flips, TR, TempRes, acq_sequence, varargin)
% Inputs:
% S - Mxy signal to fit to, size [#vox nmets Nt]
%     metabolites in order: P, L, B, A
% params_fixed - struct containing parameters (R1P, R1L, ..., Mz0P,
%     Mz0L,..., kPL, kPB, kPA) that are fixed
%     NOTE: R2s should be included for bSSFP sequence
% params_est - struct containing estimates for parameters (R1P, R1L, ..., 
%     Mz0_P, Mz0_L,..., kPL, kPB, kPA) that will be fit
%     NOTE: R2s should be included for bSSFP sequence
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
% acq_sequence - string array describing the sequence used for each
%               metabolite in orer P, L, B, A; sequence choices are:
%               "3DbSSFP", "2DGRE", "3DGRE"
%
% Optional Inputs:
% input - input function [1 Nt], if input isnt provided will estimate input
%   from pyruvate (inputless), if want no input set input to zeros([1 Nt])
% scales - vector of values to use to scale the signal after modelling of
%   size [1 Nmets] (default: [1 1 1 1])
% cat_flips - vector for catalyzation flip angle sequence for bSSFP 
%    acquisition, default is [half of 1st bSSFP flip angle]
% cat_TR - vector for catalyzation TR sequence for bSSFP acquisition,
%    default is [half of bSSFP TR] 
% B1 - vector of size [#vox 1], if B1 is provided will scale flip angles
%    with the given value per voxel
% verbose - 1 to display plots and fitting results, 0 to not, default:0
%
% Outputs:
% finalfit- a struct of parameters fit ex. finalfit.kPL corresponds to fit 
%            kPL and so on
% metrics - a struct of parameters that describe the goodness of fit
%    including objective values of lsqnonlin, NRMSE values
% Mxy_fit - an array of [#vox Nmets Nt] fit signal/ transverse mag
% Mz_fit - an array of [#vox Nmets Nt] fit longitudinal mag
% input_all - an array of [#vox Nt] for estimated input per voxel
%
% Author: Sule Sahin, 2024 Copyright

    % TO DO: add user specified lower and upper bounds  
    % TO DO: are current error metrics sufficient?     
    % TO DO: add checks in place to catch errors
    
    p = inputParser;
    p.addParameter( 'input', [], @isvector );
    p.addParameter( 'scales', [], @isvector );
    p.addParameter( 'cat_flips', [], @isvector );
    p.addParameter( 'cat_TR', [], @isvector );
    p.addParameter( 'B1', [], @isvector );
    p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
    p.parse( varargin{:} );
    input = p.Results.input;
    scales = p.Results.scales;
    cat_flips = p.Results.cat_flips;
    cat_TR = p.Results.cat_TR;
    B1vals = p.Results.B1;
    verbose = p.Results.verbose;
    
    %-------- parameter setup & variable checks --------
    if length(size(S)) == 3
        nvox = size(S,1);
        nmets = size(S,2);
        Nt = size(S,3);
    elseif length(size(S)) == 2 
        nmets = size(S,1);
        Nt = size(S,2);
        nvox = 1;
    else
        error('Size of S is incorrect-- should be [#vox nmets Nt].')
    end

    %check all variables are defined properly for sequence and set spoilers
    met_order = ["P", "L", "B", "A"];
    spoilers = ones([1 nmets]);
    for i=1:nmets
        if acq_sequence(i) == "3DbSSFP"
            spoilers(i) = 0;
            if ~any(strcmp(fieldnames(params_est), strcat("R2",met_order(i)))) && ~any(strcmp(fieldnames(params_fixed), strcat("R2",met_order(i))))
                warning(strcat(strcat("R2",met_order(i)), " not fixed or estimated but using bSSFP sequence for metabolite ",met_order(i),"."))
            end
            if size(TR.(met_order(i)),2) == 1
                error(strcat("TR for metabolite ",met_order(i)," should be size >1 as it is a 3DbSSFP acquisition."))
            end
            if size(flips.(met_order(i)),2) == 1
                error(strcat("Flips for metabolite ",met_order(i)," should be size >1 as it is a 3DbSSFP acquisition."))
            end
        elseif acq_sequence(i) == "3DGRE"
            if size(TR.(met_order(i)),2) == 1
                error(strcat("TR for metabolite ",met_order(i)," should be size >1 as it is a 3DGRE acquisition."))
            end
            if size(flips.(met_order(i)),2) == 1
                error(strcat("Flips for metabolite ",met_order(i)," should be size >1 as it is a 3DGRE acquisition."))
            end
        elseif acq_sequence(i) == "2DGRE"
            if size(TR.(met_order(i)),2) > 1
                error(strcat("TR for metabolite ",met_order(i)," should be size =1 as it is a 2DGRE acquisition."))
            end
            if size(flips.(met_order(i)),2) > 1
                error(strcat("Flips for metabolite ",met_order(i)," should be size =1 as it is a 2DGRE acquisition."))
            end
        else
            error(strcat(acq_sequence(i), " is not an accepted acquisition sequence."))
        end
    end
    
    %if no R2 specified for each metabolite set it to zero
    if ~any(strcmp(fieldnames(params_est), 'R2P')) && ~any(strcmp(fieldnames(params_fixed), 'R2P'))
        params_fixed.R2P = 0;
    end
    if ~any(strcmp(fieldnames(params_est), 'R2L')) && ~any(strcmp(fieldnames(params_fixed), 'R2L'))
        params_fixed.R2L = 0;
    end
    if ~any(strcmp(fieldnames(params_est), 'R2B')) && ~any(strcmp(fieldnames(params_fixed), 'R2B'))
        params_fixed.R2B = 0;
    end
    if ~any(strcmp(fieldnames(params_est), 'R2A')) && ~any(strcmp(fieldnames(params_fixed), 'R2A'))
        params_fixed.R2A = 0;
    end

    if isempty(input) %if no input given set as zeros and estimate per voxel later
        input_all = zeros([nvox, Nt]);
    else %if input, set same input for all voxels %%% TO DO check size of input
        input_all = repmat(input, [nvox 1]);
    end

    if isempty(scales) %if no scales given, assume no scaling in b/w mets
        scales = ones([1 nmets]);
    end
    
    if isempty(B1vals) %if no B1 provided, scale by 1
        B1vals = ones([nvox, 1]);
    else
        B1vals(B1vals==0) = 1; % to avoid scaling flips with zero
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
    pyr_cat_flips = cat_flips_tot{1};
    pyr_cat_TR = cat_TR_tot{1};

    %check TR and TempRes defined correctly
    TotalTR = 0;
    TRfields = fieldnames(TR);
    for f=1:numel(TRfields)
        TotalTR = TotalTR + sum(TR.(TRfields{f}));
    end
    TotalTR = TotalTR + sum(cat(2,cat_TR_tot{:}))*2; %add catalyzation TRs
    if abs(TempRes - TotalTR) > 0.3
        warning("Sum of TRs is not equal or close to Temporal Resolution. The difference between the sum and Temp Res is more than 0.3s. Please check definition of TR and TempRes.")
    end
    
    switch nmets
        case 2
            met_list1 = {'pyruvate', 'lactate', 'pyruvate fit', 'lactate fit'};
            met_list2 = {'pyruvate fit', 'lactate fit', 'input fit'};
            all_params = {'Mz0_P', 'Mz0_L', 'R1P', 'R1L', 'R2P', 'R2L', 'kPL'};
            all_lowb = [0, 0, 1/100, 1/100, 1/10, 1/10, 0];
            all_upb = [Inf, Inf, 1, 1, 10, 10, 1];
        case 3
            met_list1 = {'pyruvate', 'lactate', 'bicarb', 'pyruvate fit', 'lactate fit', 'bicarb fit'};
            met_list2 = {'pyruvate fit', 'lactate fit', 'bicarb fit', 'input fit'};
            all_params = {'Mz0_P', 'Mz0_L', 'Mz0_B', 'R1P', 'R1L', 'R1B', 'R2P', 'R2L', 'R2B', 'kPL', 'kPB'};
            all_lowb = [0, 0, 0, 1/100, 1/100, 1/100, 1/10, 1/10, 1/10, 0, 0];
            all_upb = [Inf, Inf, Inf, 1, 1, 1, 10, 10, 10, 1, 1];
        case 4
            met_list1 = {'pyruvate', 'lactate', 'bicarb', 'alanine', 'pyruvate fit', 'lactate fit', 'bicarb fit', 'alanine fit'};
            met_list2 = {'pyruvate fit', 'lactate fit', 'bicarb fit', 'alanine fit', 'input fit'};
            all_params = {'Mz0_P', 'Mz0_L', 'Mz0_B', 'Mz0_A', 'R1P', 'R1L', 'R1B', 'R1A', 'R2P', 'R2L', 'R2B', 'R2A', 'kPL', 'kPB', 'kPA'};
            all_lowb = [0, 0, 0, 0, 1/100, 1/100, 1/100, 1/100, 1/10, 1/10, 1/10, 1/10, 0, 0, 0];
            all_upb = [Inf, Inf, Inf, Inf, 1, 1, 1, 1, 10, 10, 10, 10, 1, 1, 1];
        otherwise
            error('Only supports 2-4 metabolites.')
    end
    
    %convert params_est to vector
    fit_params = fieldnames(params_est);
    params_est_vector = zeros([1 length(fit_params)]);
    ub = zeros([1 length(fit_params)]);
    lb = zeros([1 length(fit_params)]);
    for a=1:length(fit_params)
        if any(strcmp(all_params, fit_params{a}))
            params_est_vector(a) = params_est.(fit_params{a});
            ub(a) = all_upb(strcmp(fit_params(a),all_params));
            lb(a) = all_lowb(strcmp(fit_params(a),all_params));
        end
    end


    %-------- fitting loop---------
    fitparam = zeros([nvox, length(fit_params)]); 
    Mxy_fit = zeros([nvox, nmets, Nt]); Mz_fit = zeros([nvox, nmets, Nt]); 
    for i=1:nvox
        disp(strcat('Fitting voxel number', {' '}, string(i), ' out of ', {' '}, string(nvox)))
        if length(size(S)) == 3
            S_vox = squeeze(S(i,:,:));
        else
            S_vox = S;
        end
        
        %scale flip angles with B1 values
        mets = fieldnames(flips);
        for met=1:length(mets)
            flips_scaled.(mets{met}) = flips.(mets{met}) .* B1vals(i);
        end
        
        opts = optimoptions('lsqnonlin', 'Display', 'none', 'UseParallel', true);
        obj = @(var) difference_squared(var, params_fixed, S_vox, flips_scaled, TR, TotalTR, TempRes, pyr_cat_flips, pyr_cat_TR, cat_flips, cat_TR, Nt, spoilers, all_params, input_all(i,:), scales, fit_params) ;  % perform least-squares in signal domain
        [fitparam(i, :), metrics.obj_val] = lsqnonlin(obj, params_est_vector, lb, ub, opts);
        %[fitparam, metrics.obj_val] = lsqnonlin(obj, params_est_vector, [], [], opts);
        % TO DO: errors on kPL and other rate constants
    
        [Mz0_new, R1_new, R2_new, k_new] = unpack_parameters(fitparam(i, :), params_fixed, all_params, fit_params);
        fit_string = "|Mxy|";
        for a=1:length(fit_params)
            finalfit.(fit_params{a})(i) = fitparam(i, a);
            fit_string = strcat(fit_string, ",", " ", fit_params{a}, "=", num2str(fitparam(i, a)), " ");
        end

        if ~any(input_all(i,:)) %if input function is all zeros estimate input
            input_all(i,:) = estimate_input_function(k_new, S_vox(1,:), flips_scaled.P, R1_new(1), R2_new.P, TR.P, TotalTR, pyr_cat_flips, pyr_cat_TR, Mz0_new(1), spoilers(1), scales(1));
        end
        [Mxy_fit(i,:,:), Mz_fit(i,:,:)] = sim_multisite_bSSFP(flips_scaled, TR, TempRes, R1_new, k_new, Nt, Mz0_new, 'R2', R2_new, 'spoilers', spoilers, 'scales', scales, 'cat_flips', cat_flips, 'cat_TR', cat_TR, 'input', input_all(i,:));

        % calculate errors
        metrics.NRMSE.P(i) = rms(squeeze(Mxy_fit(i,1,:))' - S_vox(1,:))./mean(S_vox(1,:)); %pyruvate NRMSE should be ~0
        metrics.NRMSE.L(i) = rms(squeeze(Mxy_fit(i,2,:))' - S_vox(2,:))./mean(S_vox(2,:));
        if nmets>2
            metrics.NRMSE.B(i) = rms(squeeze(Mxy_fit(i,3,:))' - S_vox(3,:))./mean(S_vox(3,:));            
        end
        if nmets>3
            metrics.NRMSE.A(i) = rms(squeeze(Mxy_fit(i,4,:))' - S_vox(4,:))./mean(S_vox(4,:));
        end

        if verbose
            t=0:TempRes:(Nt-1)*TempRes;
            figure();
            subplot(211)
            plot(t, S_vox, 'LineWidth', 2)
            hold on
            plot(t, abs(squeeze(Mxy_fit(i,:,:))), '--','LineWidth', 2)
            title(fit_string)
            xlabel('time (s)')
            ylabel('|Mxy|')
            legend(met_list1)

            subplot(212)
            plot(t, squeeze(Mz_fit(i,:,:)), '--', 'LineWidth', 2)
            hold on
            plot(t, input_all, 'k--', 'LineWidth', 2)
            title('Mz')
            xlabel('time (s)')
            ylabel('Mz')
            legend(met_list2)

            %display results
            disp('----------- Fitting Results --------------')
            [~, fit_string_new] = strtok(fit_string);
            disp(fit_string_new)
            disp(strcat('Pyr NRMSE:', num2str(metrics.NRMSE.P(i))))
            disp(strcat('Lac NRMSE:', num2str(metrics.NRMSE.L(i))))
            if nmets > 2
                disp(strcat('Bicarb NRMSE:', num2str(metrics.NRMSE.B(i))))              
            end
            if nmets > 3
                disp(strcat('Ala NRMSE:', num2str(metrics.NRMSE.A(i))))
            end
        end
    end
end


function diff = difference_squared(params_est_vector, params_fixed, S, flips, TR, TotalTR, TempRes, pyr_cat_flips, pyr_cat_TR, cat_flips, cat_TR, Nt, spoilers, all_params, input, scales, fit_param_names)
    [Mz0, R1, R2, k] = unpack_parameters(params_est_vector, params_fixed, all_params, fit_param_names);
    
    if ~any(input)
        input = estimate_input_function(k, S(1,:), flips.P, R1(1), R2.P, TR.P, TotalTR, pyr_cat_flips, pyr_cat_TR, Mz0(1), spoilers(1), scales(1));
    end

    [Mxy, ~] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'R2', R2, 'spoilers', spoilers, 'scales', scales, 'cat_flips', cat_flips, 'cat_TR', cat_TR, 'input', input);
    %diff = (S(2:end,:) - Mxy(2:end,:)).^2; 
    diff = S(2:end,:) - Mxy(2:end,:); 
    %diff = (S - Mxy).^2;
    %diff = S - Mxy;
    diff = diff(:);
end


function [Mz0, R1, R2, k] = unpack_parameters(params_est_vector, params_fixed, all_params, fit_param_names)

    num_mets = sum(contains(fit_param_names, 'Mz0')) + sum(contains(fieldnames(params_fixed), 'Mz0'));
    ind1 =1; ind2 = 1; ind3=1;
    Mz0 = zeros(num_mets, 1); k=zeros(1, num_mets-1); R1 = zeros(num_mets, 1);
    
    for a=1:length(all_params)
        if any(strcmp(fit_param_names, all_params{a}))
            if contains(all_params{a}, 'Mz0')
                Mz0(ind1) = params_est_vector(strcmp(fit_param_names, all_params{a}));
                ind1 = ind1+1;
            elseif contains(all_params{a}, 'R1') 
                R1(ind2) = params_est_vector(strcmp(fit_param_names, all_params{a}));
                ind2 = ind2+1;
            elseif contains(all_params{a}, 'R2') 
                met = all_params{a}(end);
                R2.(met) = params_est_vector(strcmp(fit_param_names, all_params{a}));
            elseif contains(all_params{a}, 'kP') 
                k(ind3) = params_est_vector(strcmp(fit_param_names, all_params{a}));
                ind3 = ind3+1;
            end
        elseif isfield(params_fixed, all_params{a})
            if contains(all_params{a}, 'Mz0')
                Mz0(ind1) = params_fixed.(all_params{a});
                ind1 = ind1+1;
            elseif contains(all_params{a}, 'R1') 
                R1(ind2) = params_fixed.(all_params{a});
                ind2 = ind2+1;
            elseif contains(all_params{a}, 'R2') 
                met = all_params{a}(end);
                R2.(met) = params_fixed.(all_params{a});
            elseif contains(all_params{a}, 'kP') 
                k(ind3) = params_fixed.(all_params{a});
                ind3 = ind3+1;
            end
        else
            error(strcat('Did not include value for parameter below',str(all_params{a})))
        end
    end
end


function input_fxn = estimate_input_function(k, Pxy, pyrflips, R1P, R2P, TRP, TotalTR, pyr_cat_flips, pyr_cat_TR, Mz0P, spoilers, norm)
    
    TempRes_nopyr = TotalTR - (sum(TRP)+2*sum(pyr_cat_flips));
    
    if size(pyrflips, 2) == 1
        input_fxn = input_est_2D(pyrflips, R1P, R2P, k, TRP, TempRes_nopyr, Mz0P, spoilers, Pxy, norm);
    else
        input_fxn = input_est_3D(pyrflips, R1P, R2P, k, TRP, TempRes_nopyr, pyr_cat_flips, pyr_cat_TR, Mz0P, spoilers, Pxy, norm);
    end
end

function [input] = input_est_2D(flips, R1, R2, k, TR, TempRes, Mz0P, spoilers, Pxy, norm)
    Nt = size(Pxy,2);
    Pz = (Pxy ./ norm) ./ sind(flips)'; 
    spoil = diag([0 0 1])^spoilers;
    R = get_R(R1, R2, k);
    
    input = zeros([1, Nt]);
    for nt=1:Nt-1
        M_t = rot_matrix(flips(nt,:)) * spoil * [0;0;Pz(nt)];
        M_next = spoil * expm(R*(TR+TempRes)) * M_t; %mag before next RF
        input(nt+1) = Pz(nt+1) - M_next(3);
    end
    input(1) = Pz(1) - Mz0P;
end

function [input] = input_est_3D(flips_nocat, R1, R2, k, TR_nocat, TempRes, cat_flips, cat_TR, Mz0P, spoilers, Pxy, norm)
    
    Nt = size(Pxy,2);
    spoil = diag([0 0 1])^spoilers;
    R = get_R(R1, R2, k);

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
    
    % convert Pxy to Pz
    scale = zeros([1 Nt]);
    for s=1:Nt
        M_ti = zeros([3 3 size(flips_nocat,2)]);
        for i=1:num_ex
            rot = rot_matrix(flips(s,i));
            if i==1
                M_next = rot;
            else
                M_next = rot * expm(R*TR(i-1)) * spoil * M_next;
            end
            if and(i >= start_idx, i <= end_idx)
                M_ti(:,:,i-cat_len) = M_next;
            end
        end
        scale(s) = abs(mean(M_ti(1,3,:),3));
    end
    
    Pz = (Pxy ./ norm) ./ scale; %assume each time point is the same
    %%% TO DO: make this work with variation per time point
    
    input = zeros([1, Nt]);
    for nt=1:Nt-1
        for i=1:num_ex
            rot = rot_matrix(flips(nt, i));
            if i==1
                M_next = rot * spoil * [0;0;Pz(nt)];
            else
                M_next = rot * expm(R*TR(i-1)) * spoil * M_next;
            end
        end
        M_next = expm(R*(TR(end)+TempRes)) * spoil * M_next; %last bit of relaxation after last flip
        input(nt+1) = Pz(nt+1) - M_next(3);
    end
    input(1) = Pz(1) - Mz0P;
    
end


function R = get_R(R1, R2, k)
    R = [-R2  0      0      ;
           0  -R2    0      ;
           0    0  -sum(k)-R1];
end

function rot = rot_matrix(FA)

    rot = [ cosd(FA) 0 -sind(FA);
        0 1 0;
        sind(FA) 0 cosd(FA)];

end