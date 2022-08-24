close all; clear all;

cd /Users/sf565268/Documents/matlab/SSFP/L+S_simulations/

for aa = 1:1

load epi_TRAMP_tumor.mat

for i = 1:16
    tumor1(:,:,i) = pa(:,:,i,8);
    tumor(:,:,i) = imresize(tumor1(3:12,3:16,i),[20,20]);
end

data = ifftnc_time(tumor,3);
data = data/max(data(:));
% data = imnoise(data,'gaussian',0,0.0000000005);
pa = fftnc_time(data,3);
img_noise = std2(abs(pa(1:5,1:5,1)));
img_SNR = 0.655*pa/img_noise;
y = size(data);
% [u,v,t] = svd(reshape(epi_kidney_reformat,[y(1)*y(2) y(3)]),0);
% data = permute(data,[2 3 1]);
%% Undersampling Pattern
Sampling_Percentage = 25;
str = ['pattern_' int2str(Sampling_Percentage) '_percent_tramp.mat'];
load (str);


undersample = pattern;
% undersample = ones(y(1),y(2),y(3));

kdata = data.*undersample;

% kdata = permute(kdata,[3 1 2]);
% data = permute(data,[3 1 2]);
undersample = zeros(y(1),y(2),y(3));
for i=1:y(1)
    for j=1:y(2)
        
            for l=1:y(3)
                if abs(kdata(i,j,l))>0
                    undersample(i,j,l) = 1;
                elseif abs(kdata(i,j,l)) == 0
                    undersample(i,j,l)=0;
                end
            end
        
    end
end
    
DN = [y(1) y(2) y(3)];
timedim = 3; 
wavelet_time_scale = 1;
sparse = 0;
nIter = 1;
LplusS = zeros(y(1),y(2),y(3),nIter);
%% L+S reconstruction **************
tic
if sparse == 0;
    for i=1:nIter
    param.E=(pnDFT_time(undersample,DN,timedim,1,2));
    param.d=kdata;
%     param.T=Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);
    param.T= TempFFT(3);
    param.lambda_L=0.01;
    param.lambda_S=0.001;
    param.nite=800;
    param.tol=1.5e-3;
    param.bsize = [5,5];
    param.mu = 2;
    param.hard_thresh = 0.3;

    [L,S] = LLR_plus_S_PCA(param); 

    LplusS(:,:,:,i)=L+S;
    end
    
else
    param.E=pnDFT_time(undersample,DN,timedim,1,2);
    param.W = TempFFT(3);
%     param.W=Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);
    param.L1Weight=0.0004;
    param.TV = TVOP_3D();param.TVWeight=0.0000008;
    param.y = kdata;
    param.nite = 24;
    param.display=1;

    % linear reconstruction
    recon_dft=param.E'*kdata;

    % k-t SPARSE-SENSE reconstruction
    fprintf('\n k-t SPARSE-SENSE reconstruction \n')
    tic
    recon_cs=recon_dft;
    % the non-linear CG algorithm will re-start two times to improve
    % convergence using the last result as starting point
    for h=1:4,
        recon_cs = CSL1NlCg(recon_cs,param); L=zeros(y(1),y(2),y(3)); S=zeros(y(1),y(2),y(3));
%         recon_cs = CSL1NlCg_pca(recon_cs,param); L=zeros(y(1),y(2),y(3)); S=zeros(y(1),y(2),y(3));
    end
    LplusS = recon_cs;
    toc
end
toc
%%
LplusS = mean(LplusS,4);
pa = fftnc_time(data,3);
LplusS = LplusS/max(LplusS(:));
pa = pa/max(pa(:));
recon_dft = param.E'*kdata;
recon_dft = recon_dft/max(recon_dft(:));
for i = 1:y(3)
    LplusS(:,:,i) = flipud(LplusS(:,:,i));
%     L(:,:,i) = imtranslate(flipud(L(:,:,i)),[1 -1]);
%     S(:,:,i) = imtranslate(flipud(S(:,:,i)),[1 -1]);
%     LplusS(:,:,i) = imtranslate(LplusS(:,:,i),[1 -1]);
    LplusS(:,:,i) = circshift(LplusS(:,:,i),[1 -1]);
    pa(:,:,i) = fliplr(pa(:,:,i));
    recon_dft(:,:,i) = flipud(recon_dft(:,:,i));
    recon_dft(:,:,i) = circshift(recon_dft(:,:,i),[1 -1]);
%     diff(:,:,i) = imtranslate(diff(:,:,i),[-1 1]);
end
differ = abs(LplusS)-abs(pa);

%%
figure;imagesc(abs(pa(:,:,4))); axis image; set(gca,'xtick',[],'ytick',[])
figure;imagesc(abs(LplusS(:,:,4))); axis image; set(gca,'xtick',[],'ytick',[])
figure;imagesc(abs(recon_dft(:,:,4))); axis image; set(gca,'xtick',[],'ytick',[])
figure;imagesc(abs(differ(:,:,4)),[0 1]); axis image; set(gca,'xtick',[],'ytick',[])

nrmse = RMSE(abs(LplusS(:)),abs(pa(:)))/(max(abs(pa(:)))-min(abs(pa(:))));
nrmse2 = RMSE(abs(LplusS(:)),abs(pa(:)))/(rms(abs(pa(:))));
nrmse3 = RMSE(LplusS(:),pa(:))/(rms(pa(:)));

%%
ssimval = zeros(1,y(3));
ssim_map = zeros(y(1),y(2),y(3));
for i=1:y(3)
    [ssimval(1,i),ssim_map(:,:,i)] = ssim_index(abs(LplusS(:,:,i)),abs(pa(:,:,i)));
end
mean_ssimval = mean(ssimval(1,3:end));
nRMSE = norm(LplusS(:) - pa(:)) / norm(pa(:));
nRMSE_abs = norm(abs(LplusS(:)) - abs(pa(:))) / norm(abs(pa(:)));
nRMSE_abs_dft = norm(abs(recon_dft(:)) - abs(pa(:))) / norm(abs(pa(:)));
%%
vasculature_LplusS = zeros(1,y(3));
vasculature_orginal = zeros(1,y(3));
for i=1:y(3)
    vasculature_LplusS(1,i) = mean2(abs(LplusS(8:16,8:12,i)));
    vasculature_original(1,i) = mean2(abs(pa(8:16,8:12,i)));
    vasculature_dft(1,i) = mean2(abs(recon_dft(8:16,8:12,i)));
end
x = 1:y(3);
figure;
plot(x,vasculature_original,x,vasculature_LplusS,'--',x,vasculature_dft,'LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('Ground Truth','LplusS','Zero-filled Recon','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')
%%
kidney_LplusS = zeros(1,y(3));
kidney_orginal = zeros(1,y(3));
for i=1:y(3)
    kidney_LplusS(1,i) = mean2(abs(LplusS(4:13,3:10,i)));
    kidney_original(1,i) = mean2(abs(pa(4:13,3:10,i)));
    kidney_dft(1,i) = mean2(abs(recon_dft(4:13,3:10,i)));
end
x = 1:y(3);
figure;
plot(x,kidney_original,x,kidney_LplusS,'--',x,kidney_dft,'LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('Ground Truth','LplusS','Zero-filled Recon','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')

%%
[u_pa,v_pa,t_pa] = svd(reshape(pa,[y(1)*y(2) y(3)]),0);
v_pa = diag(v_pa)/max(v_pa(:));
[u_lps,v_lps,t_lps] = svd(reshape(LplusS,[y(1)*y(2) y(3)]),0);
v_lps = diag(v_lps)/max(v_lps(:));
[u_dft,v_dft,t_dft] = svd(reshape(recon_dft,[y(1)*y(2) y(3)]),0);
v_dft = diag(v_dft)/max(v_dft(:));
x = 1:y(3);
figure;
semilogy(x,v_pa,x,v_lps,x,v_dft)
%{
%%
figure;
for i=1:y(3)
    raw_image = imresize(abs(LplusS(:,:,i)), [64 64]);
    raw_image = rot90(raw_image);
    subplot_tight(1,y(3),i,0.001); imagesc(raw_image,[0 1]); colormap(jet);  axis off;  
end

%%
figure;
for i=1:y(3)
    raw_image = imresize(abs(pa(:,:,i)), [64 64]);
    raw_image = rot90(raw_image);
    subplot_tight(1,y(3),i,0.001); imagesc(raw_image,[0 1]); colormap(jet);  axis off;  
end

%%
figure;
for i=1:y(3)
    raw_image = imresize(ssim_map(:,:,i), [64 64]);;
    raw_image = rot90(raw_image);
    subplot_tight(1,y(3),i,0.001); imagesc(raw_image,[0 1]); colormap(jet);  axis off;  
end
%}
%%
LplusS_notnorm = L+S;
L = L/max(L(:));
S = S/max(S(:));
for i=1:y(3)
    vasculature_LplusS_notnorm(1,i) = mean2(abs(LplusS_notnorm(15:19,7:10,i)));
    vasculature_L(1,i) = mean2(abs(L(15:19,7:10,i)));
    vasculature_S(1,i) = mean2(abs(S(15:19,7:10,i)));
end
x = 1:y(3);
figure;
plot(x,vasculature_LplusS_notnorm,x,vasculature_L,'--',x,vasculature_S,'LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('LplusS','L','S','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')

for i=1:y(3)
    kidney_LplusS_notnorm(1,i) = mean2(abs(LplusS_notnorm(4:13,3:10,i)));
    kidney_L(1,i) = mean2(abs(L(4:13,3:10,i)));
    kidney_S(1,i) = mean2(abs(S(4:13,3:10,i)));
end
x = 1:y(3);
figure;
plot(x,kidney_LplusS_notnorm,x,kidney_L,'--',x,kidney_S,'LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('LplusS','L','S','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')

if aa <10
	close all;
end

ssimval2(1,aa) = mean_ssimval;
nRMSE2(1,aa) = nrmse2;
aa
end
mean_ssimval2 = mean(ssimval2);
mean_nRMSE2 = mean(nRMSE2);