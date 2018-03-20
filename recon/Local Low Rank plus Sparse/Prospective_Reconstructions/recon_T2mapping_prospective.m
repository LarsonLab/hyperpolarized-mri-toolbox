clear all;
clc;
close all;

load t2mapping_c2pyr_prospective_data
load t2mapping_lactate_prospective_data
load t2mapping_HP001_prospective_data
load t2mapping_urea_prospective_data

data = squeeze(data);
% data = data(:,:,:);
data = data/max(data(:));
y = size(data);

mask = abs(data)>0;
DN = [y(1) y(2) y(3)];
timedim = 3; 
wavelet_time_scale = 1;
sparse = 0;
nIter = 1;
LplusS = zeros(y(1),y(2),y(3),nIter);
%%
tic
if sparse == 0;
    for i=1:nIter
    param.E=pnDFT_time(mask,DN,timedim,1,2);
    param.d=data;
%     param.T=Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);
    param.T= TempFFT(3);
    param.lambda_L=0.1;
    param.lambda_S=0.001;
    param.nite=500;
    param.tol= 4e-3;
%     param.tol=2.3e-3;
    param.bsize = [5,5];
    param.mu = 2;
    param.hard_thresh = 0.3;

    [L,S] = LLR_plus_S_PCA(param); 

    LplusS(:,:,:,i)=L+S;
    end
else
    param.E=pnDFT_time(mask,DN,timedim,1,2);
%     param.W = TempFFT(3);
%     param.W = coeff;
%     param.W=Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);
    param.L1Weight=0.01;
    param.TV = TVOP_3D();param.TVWeight=0.000000;
    param.y = data;
    param.nite = 30;
    param.display=1;

    % linear reconstruction
    recon_dft=param.E'*data;

    % k-t SPARSE-SENSE reconstruction
    fprintf('\n k-t SPARSE-SENSE reconstruction \n')
    tic
    recon_cs=recon_dft;
    % the non-linear CG algorithm will re-start two times to improve
    % convergence using the last result as starting point
    for h=1:10,
        recon_cs = CSL1NlCg_pca(recon_cs,param);
    end
    LplusS = recon_cs;
    toc
end
toc

%%
LplusS = mean(LplusS,4);
for i = 1:size(LplusS,3)
    LplusS(:,:,i) = flipud(fliplr(LplusS(:,:,i)));
end
%%
figure;imagesc(abs(LplusS(:,:,1)))
kspace = ifftnc_time(LplusS,3);
figure;imagesc(abs(kspace(:,:,1)))

%%
recon_dft = param.E'*data;
LplusS = LplusS/max(LplusS(:));
recon_dft = recon_dft/max(recon_dft(:));
diff = abs(LplusS)-abs(recon_dft);