function [data_Out,varargout] = pb_phase_corr(data_In,spec_idx_met,varargin)
% convert spectra to image, sum over met peaks
% v0: 20180816
    calc_phi_0_flag = 1;
    if nargin > 2
        calc_phi_0_flag = 0;
        phi_0 = varargin{1};
        fprintf('Correct zero-order phase based on table ... \n');
    end
    size_Data_5D = size(data_In);
    
    % apply linear phase
    phi_1 = 3.92*pi; % manually tuned for the sequence
    f = linspace(-1/2,1/2,size_Data_5D(1));
    data_In = data_In .* exp(-1i*phi_1*repmat(f',[1 size_Data_5D(2:4)]));
    
    % detect 0th order phase for each resonance
    if calc_phi_0_flag
        Data_5D_cplxAUC = sum(data_In,4);
        for i_x = 1:size_Data_5D(2),
        for i_y = 1:size_Data_5D(3),
        for i_mets = 1:size(spec_idx_met,2),
            phi_0(i_x,i_y,i_mets) = find_phase_corr(...
                Data_5D_cplxAUC(spec_idx_met{i_mets},i_x,i_y));
        end
        end
        end
    end
    % zero-order phase correction
    data_Out = data_In .* exp(1i*repmat(shiftdim(phi_0(:,:,1),-1),...
        [size_Data_5D(1) 1 1 size_Data_5D(4)]));
    varargout{1} = phi_0;
end