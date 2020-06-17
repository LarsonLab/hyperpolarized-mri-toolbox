function dataOut = pb_splicing(dataIn,dataIn_bicarb,spec_idx_met)

    size_dataIn = size(dataIn);
    length_lac = length(spec_idx_met{2});
    length_bicarb = length(spec_idx_met{3});
    dataOut = zeros([length_bicarb+length_lac size_dataIn(2:end)]);
    dataOut(1:length_bicarb,:,:,:) = dataIn_bicarb(spec_idx_met{3},:,:,:);
    dataOut((length_bicarb+1):end,:,:,:) = dataIn(spec_idx_met{2},:,:,:);

end