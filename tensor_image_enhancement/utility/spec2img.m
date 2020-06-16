function met_specAUC = spec2img(Data_5D,spec_idx_met)
% convert spectra to image, sum over met peaks
% v0: 20180816
    % detect 0th order phase for each resonance
    size_Data_5D = size(Data_5D);
    Data_5D_cplxAUC = sum(Data_5D,4);
    for i_x = 1:size_Data_5D(2),
    for i_y = 1:size_Data_5D(3),
    for i_mets = 1:size(spec_idx_met,2),
        phi_0(i_x,i_y,i_mets) = find_phase_corr(...
            Data_5D_cplxAUC(spec_idx_met{i_mets},i_x,i_y));
    end
    end
    end

    % phase correction and peak estimation
    for i_mets = 1:size(spec_idx_met,2),
    met_spec = Data_5D .* exp(1i*repmat(shiftdim(phi_0(:,:,i_mets),-1),...
        [size_Data_5D(1) 1 1 size_Data_5D(4)]));
    met_specAUC{i_mets} = squeeze(sum(real(met_spec(spec_idx_met{i_mets},:,:,:))));
    end

end