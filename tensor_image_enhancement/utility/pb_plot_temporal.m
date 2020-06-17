function pb_plot_temporal(met_specAUC_before,met_specAUC_after)
% v0: plot time series 2DEPSI
% v2: use montage
    name_mets = {'Pyruvate','Lactate','Bicarbonate'};
    for i_mets = 1:size(met_specAUC_before,2)
        met_specAUC{i_mets} = cat(3,met_specAUC_before{i_mets},met_specAUC_after{i_mets});
        met_cmax(i_mets) = max(met_specAUC{i_mets}(:));
    end
%     met_cmax = [1.8e9 2.1e8];
    figure,
    for i_mets = 1:size(met_specAUC_before,2)
        subplot(3,1,i_mets)
        temp1 = permute(flip(met_specAUC{i_mets},1),[2,1,4,3]);
        montage(temp1, 'Size', [2 size(met_specAUC_before{i_mets},3)], 'DisplayRange', [0 met_cmax(i_mets)]);
        axis image off, colormap(gray(256));
        title(name_mets{i_mets},'Fontsize',14);
    end

%     set(subplot(2,1,1), 'Position', [0.05, 0.69, 0.92, 0.27])
%     set(subplot(2,1,2), 'Position', [0.05, 0.37, 0.92, 0.27])
    set(subplot(3,1,1), 'Position', [0.05, 0.69, 0.92, 0.27])
    set(subplot(3,1,2), 'Position', [0.05, 0.37, 0.92, 0.27])
    set(subplot(3,1,3), 'Position', [0.05, 0.05, 0.92, 0.27])

%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.7, 0.8]);


end