function [CF_Params,TdCostFxn] = td_CostFxnTrade(TdBias,TdVariance,BVratio)
% tradeoff bias and variance through cost function v1. 20190207
    TdCostFxn = TdBias + BVratio*TdVariance;
    [y_CostFxn i_CostFxn] = sort(TdCostFxn(:),'ascend');
    % pick out the terms where bias > bias @ mincost (stronger denoising)
    TdBias_temp = TdBias(i_CostFxn);
    TdBias_large_logic = (TdBias_temp > TdBias_temp(1));
    TdBias_small_logic = (TdBias_temp <= TdBias_temp(1));
    i_BiasPriority = [flip(i_CostFxn(TdBias_large_logic),1);i_CostFxn(TdBias_small_logic)];
    [CFSpec,CFX,CFY,CFDyn] = ind2sub(size(TdCostFxn),i_BiasPriority);
    CF_Params.SortedRank = [CFSpec CFX CFY CFDyn];
    CF_Params.IdxMinCost = sum(TdBias_large_logic)+1;


end