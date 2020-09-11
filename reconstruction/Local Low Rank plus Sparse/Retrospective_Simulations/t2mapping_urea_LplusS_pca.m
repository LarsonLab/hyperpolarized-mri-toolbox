close all; 
clear all;

load t2mapping_urea_retrospective_data

data=squeeze(data(:,:,:));
images = fftnc_time(data,3);

data = data/max(data(:));
y = size(data);
% [u,v,t] = svd(reshape(epi_kidney_reformat,[y(1)*y(2) y(3)]),0);
data = permute(data,[2 3 1]);
%% Undersampling Pattern

load pattern_50_20percent_70x20

undersample = pattern(:,:);
% undersample = ones(y(2),y(3));
for i=1:y(1)
    kdata(:,:,i) = data(:,:,i).*undersample;
end
data = permute(data,[3 1 2]);
kdata = permute(kdata,[3 1 2]);
% undersample = pattern;

% [coeff] = pca(reshape(kdata, [y(1)*y(2) y(3)]));

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
    param.lambda_L=0.1;
    param.lambda_S=0.001;
    param.nite=500;
    param.tol=4e-3;
%     param.tol=2.3e-3;
    param.bsize = [5,5];
    param.mu = 2;
    param.hard_thresh = 0.3;

    [L,S] = LLR_plus_S_PCA(param); 

    LplusS(:,:,:,i)=L+S;
    
    end
    
else
    param.E=pnDFT_time(undersample,DN,timedim,1,2);
%     param.W = TempFFT(3);
%     param.W = coeff;
%     param.W=Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);
    param.L1Weight=0.001;
    param.TV = TVOP_3D();param.TVWeight=0.0000008;
    param.y = kdata;
    param.nite = 30;
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
        recon_cs = CSL1NlCg_pca(recon_cs,param);
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
figure;imagesc(differ(:,:,4),[0 1]); axis image; set(gca,'xtick',[],'ytick',[])
nrmse = RMSE(abs(LplusS(:)),abs(pa(:)))/(max(abs(pa(:)))-min(abs(pa(:))));
nrmse2 = RMSE(abs(LplusS(:)),abs(pa(:)))/(rms(abs(pa(:))));
nrmse3 = RMSE(LplusS(:),pa(:))/(rms(pa(:)));

%%
ssimval = zeros(1,y(3));
ssim_map = zeros(34,49,y(3));
for i=1:y(3)
    [ssim_map(:,:,i)] = ssim_index(abs(LplusS(34:67,13:61,i)),abs(pa(34:67,13:61,i)));
    ssimval(1,i) = mean2(ssim_map(:,:,i));
end
mean_ssimval = mean(ssimval(1,1:12))
nRMSE = norm(LplusS(:) - pa(:)) / norm(pa(:));
nRMSE_abs = norm(abs(LplusS(:)) - abs(pa(:))) / norm(abs(pa(:)));
nRMSE_abs_dft = norm(abs(recon_dft(:)) - abs(pa(:))) / norm(abs(pa(:)));
% LplusS = mean(LplusS,4);
% pa = images;
% LplusS = LplusS/max(LplusS(:));
% pa = pa/max(pa(:));
% recon_dft = param.E'*kdata;
% recon_dft = recon_dft/max(recon_dft(:));
% diff = LplusS-pa;
% for i = 1:y(3)
%     LplusS(:,:,i) = flipud(LplusS(:,:,i));
% %     L(:,:,i) = imtranslate(flipud(L(:,:,i)),[1 -1]);
% %     S(:,:,i) = imtranslate(flipud(S(:,:,i)),[1 -1]);
% %     LplusS(:,:,i) = imtranslate(LplusS(:,:,i),[1 -1]);
%     LplusS(:,:,i) = circshift(LplusS(:,:,i),[1 -1]);
%     pa(:,:,i) = fliplr(pa(:,:,i));
%     recon_dft(:,:,i) = flipud(recon_dft(:,:,i));
%     recon_dft(:,:,i) = circshift(recon_dft(:,:,i),[1 -1]);
% %     diff(:,:,i) = imtranslate(diff(:,:,i),[-1 1]);
% end
% 
% 
% % display 4 frames
% recon_cs = LplusS;
% images3D = pa;
% 
% figure;imagesc(abs(pa(:,:,4))); axis image
% figure;imagesc(abs(LplusS(:,:,4))); axis image
% figure;imagesc(abs(recon_dft(:,:,4))); axis image
% % recon_dft2=recon_dft(:,:,3);recon_dft2=cat(2,recon_dft2,recon_dft(:,:,7));recon_dft2=cat(2,recon_dft2,recon_dft(:,:,10));recon_dft2=cat(2,recon_dft2,recon_dft(:,:,y(3)));
% % recon_cs2=recon_cs(:,:,3);recon_cs2=cat(2,recon_cs2,recon_cs(:,:,7));recon_cs2=cat(2,recon_cs2,recon_cs(:,:,10));recon_cs2=cat(2,recon_cs2,recon_cs(:,:,y(3)));
% % images3D2=images3D(:,:,3);images3D2=cat(2,images3D2,images3D(:,:,7));images3D2=cat(2,images3D2,images3D(:,:,10));images3D2=cat(2,images3D2,images3D(:,:,y(3)));
% 
% 
% recon_dft2=recon_dft(:,:,3);recon_dft2=cat(2,recon_dft2,recon_dft(:,:,7));
% recon_cs2=recon_cs(:,:,3);recon_cs2=cat(2,recon_cs2,recon_cs(:,:,7));
% images3D2=images3D(:,:,3);images3D2=cat(2,images3D2,images3D(:,:,7));
% figure;
% subplot(3,1,1),imshow(abs(recon_dft2),[]);title('Zero-filled FFT')
% subplot(3,1,2),imshow(abs(recon_cs2),[]);title('Recon')
% subplot(3,1,3),imshow(abs(images3D2),[]);title('Original')
% colormap(jet)
% 
% diff = recon_cs-images3D;
% nrmse = RMSE(abs(LplusS(:)),abs(pa(:)))/(max(abs(pa(:)))-min(abs(pa(:))));
% nrmse2 = RMSE(abs(LplusS(:)),abs(pa(:)))/(rms(abs(pa(:))));
% rmse_dft = RMSE(abs(recon_dft(:)),abs(pa(:)))/(max(abs(pa(:)))-min(abs(pa(:))));
% NMSE = ((norm((LplusS(:)-pa(:)),2))^2)/(norm(pa(:),2)^2)
% 
% %%
% ssimval = zeros(1,y(3));
% ssim_map = zeros(y(1),y(2),y(3));
% for i=1:y(3)
%     [ssimval(1,i),ssim_map(:,:,i)] = ssim_index(abs(LplusS(:,:,i)),abs(pa(:,:,i)));
% end
% mean_ssimval = mean(ssimval(1,3:end));
% nRMSE = norm(LplusS(:) - pa(:)) / norm(pa(:));
% nRMSE_abs = norm(abs(LplusS(:)) - abs(pa(:))) / norm(abs(pa(:)));
% nRMSE_abs_dft = norm(abs(recon_dft(:)) - abs(pa(:))) / norm(abs(pa(:)));
%%
vasculature_LplusS = zeros(1,y(3));
vasculature_orginal = zeros(1,y(3));
for i=1:y(3)
    vasculature_LplusS(1,i) = mean2(abs(LplusS(45:48,21:24,i)));
    vasculature_original(1,i) = mean2(abs(pa(45:48,21:24,i)));
    vasculature_dft(1,i) = mean2(abs(recon_dft(45:48,21:24,i)));
end
x = 1:y(3);
figure;
plot(x,vasculature_original,x,vasculature_LplusS,'--','LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('Ground Truth','LplusS','Zero-filled Recon','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')
%%
kidney_LplusS = zeros(1,y(3));
kidney_orginal = zeros(1,y(3));
for i=1:y(3)
    kidney_LplusS(1,i) = mean2(abs(LplusS(60:64,48:52,i)));
    kidney_original(1,i) = mean2(abs(pa(60:64,48:52,i)));
    kidney_dft(1,i) = mean2(abs(recon_dft(60:64,48:52,i)));
end
nRMSE_curves = norm(abs(kidney_LplusS(:))-abs(kidney_original(:)))/norm(abs(kidney_original(:)))
x = 1:y(3);
figure;
plot(x,kidney_original,x,kidney_LplusS,'--','LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('Ground Truth','LLR+S','Zero-filled Recon','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')
%{
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
%}

%% T2mapping
infile = ifftnc_time(images,3);

% infile = kspace(:,:,1:20);

nsvd = 7;

tr = 0.013;
snrcutoff = 40;
nskip = 0;
t2lims = [0.3 20];
nt2 = 128;
lambda = 1e-4;
timepoints = size(infile,3);
[t2,t2spectra,imgs,t] = t2_spectra_mapping_ssfp_compsense_V2_lambdamax(infile,tr,snrcutoff,nskip,t2lims,nt2,timepoints,nsvd);

y = size(imgs);
t2simgs = zeros(y(1),y(2));
t2samps = zeros(y(1),y(2));
t2limgs = zeros(y(1),y(2));
t2lamps = zeros(y(1),y(2));
t2meanimgs = zeros(y(1),y(2));
t2meanamps = zeros(y(1),y(2));


for i = 1:y(1)
    for j = 1:y(2)
        for k = 1:nt2
            if t2spectra(i,j,k) > 0 && t2(1,k) <= 2.02
%                 [m,I] = mean(t2spectra(i,j,4:73));
%                 t2samps(i,j) = sum(t2spectra(i,j,4:73))/nnz(t2spectra(i,j,4:73));
                t2samps(i,j) = sum(t2spectra(i,j,1:58));
                t2simgs(i,j) = exp(sum(log(t2(1,1:58)).*squeeze(t2spectra(i,j,1:58))')./sum(squeeze(t2spectra(i,j,1:58))'));
            else
                continue
            end
        end
    end
end
for i = 1:y(1)
    for j = 1:y(2)
        for k = 1:nt2
            if t2spectra(i,j,k) > 0 && t2(1,k) > 2.03 && t2(1,k) < 41
%                 t2lamps(i,j) = sum(t2spectra(i,j,74:121))/nnz(t2spectra(i,j,74:121));
                t2lamps(i,j) = sum(t2spectra(i,j,59:nt2));
                t2limgs(i,j) = exp(sum(log(t2(1,59:nt2)).*squeeze(t2spectra(i,j,59:nt2))')./sum(squeeze(t2spectra(i,j,59:nt2))'));
            else
                continue
            end
        end
    end
end
for i = 1:y(1)
    for j = 1:y(2)
        t2meanamps(i,j) = sum(t2spectra(i,j,1:128));
        t2meanimgs(i,j) = exp(sum(log(t2(1,1:128)).*squeeze(t2spectra(i,j,1:128))')./sum(squeeze(t2spectra(i,j,1:128))'));
    end
end
                                  
t2meanimgs(isnan(t2meanimgs))=0;
close all;
%% Plotting individual short and long T2 components and amplitudes

fs = 15;
% Short T2 Component Image
h = figure;
imagesc(t2simgs,[0 2])
set(gca,'XTick',[],'YTick',[])
set(gca,'fontsize',fs)
j = jet;
j(1,:) = [0 0 0];
colorbar;
colormap(j);
ylabel(colorbar,'T_{2S}(s)','fontsize',20);
% title('[1-^{13}C]Lactate')
title('[1-13C]Pyruvate')
% saveas(h,'Data_20150724_pyruvate_t2_short','fig');

% Long T2 Component Image
i=figure;
imagesc(t2limgs,[0 10])
set(gca,'XTick',[],'YTick',[])
set(gca,'fontsize',fs)
j = jet;
j(1,:) = [0 0 0];
colorbar;
colormap(j);
ylabel(colorbar,'T_{2L}(s)','fontsize',20);
title('[1-^{13}C]Lactate')
% title('[1-13C]Pyruvate')
% saveas(i,'Data_20150724_pyruvate_t2_long','fig');

fs = 15;
% Short T2 Component Image Amplitudes
m=figure;
imagesc(t2samps,[0 100])
set(gca,'XTick',[],'YTick',[])
set(gca,'fontsize',fs)
j = jet;
j(1,:) = [0 0 0];
colorbar;
colormap(j);
ylabel(colorbar,'Amplitude','fontsize',20);
title('[1-^{13}C]Lactate')
% title('[1-13C]Pyruvate')
% saveas(m,'Data_20150724_pyruvate_t2_short_amps','fig');

% Long T2 Component Image Amplitudes
k=figure;
imagesc(t2lamps,[0 100])
set(gca,'XTick',[],'YTick',[])
set(gca,'fontsize',fs)
j = jet;
j(1,:) = [0 0 0];
colorbar;
colormap(j);
ylabel(colorbar,'Amplitude','fontsize',20);
title('[1-^{13}C]Lactate')
% title('[1-13C]Pyruvate')
% saveas(k,'Data_20150724_pyruvate_t2_long_amps','fig');

%% Plotting logarithmic mean
z=figure;
imagesc(t2meanimgs,[0 5])
set(gca,'XTick',[],'YTick',[])
set(gca,'fontsize',15)
j = jet;
j(1,:) = [0 0 0];
colorbar;
colormap(j);
ylabel(colorbar,'T_{2} mean (s)','fontsize',20);
% saveas(z,'Data_20150805_lactate_t2_mean','fig');
%% Looking at R2 Value
r2 = zeros(y(1),y(2));
for ii = 1:y(1)
    for jj = 1:y(2)
        if t2limgs(ii,jj)==0
            t2limgs(ii,jj)=1;
        end
        if t2simgs(ii,jj)==0
            t2simgs(ii,jj)=1;
        end  
        modelfit = t2lamps(ii,jj).*exp(-t/t2limgs(ii,jj))+t2samps(ii,jj).*exp(-t/t2simgs(ii,jj));
        decaycurve = squeeze(imgs(ii,jj,:));
        r2(ii,jj) = rsquare(decaycurve',modelfit);
    end
end
% for ii = 1:y(1)
%     for jj = 1:y(2)
%         modelfit = t2meanamps(ii,jj)*exp(-t/t2meanimgs(ii,jj));
%         decaycurve = squeeze(imgs(ii,jj,:));
%         r2(ii,jj) = rsquare(decaycurve',modelfit);
%     end
% end

figure;
imagesc(r2)
sum(sum(r2))/nnz(r2)

for ii = 1:y(1)
    for jj = 1:y(2)
        if r2(ii,jj) > 0.9
            continue
        else
            r2(ii,jj) = 0;
        end
    end
end
r2_mask = zeros(y(1),y(2));
for ii = 1:y(1)
    for jj = 1:y(2)
        if r2(ii,jj) > 0
            r2_mask(ii,jj) = 1;
        else
            r2_mask(ii,jj) = 0;
        end
    end
end
figure;
imagesc(r2)
figure;
imagesc(r2_mask)

t2_mean_masked = r2_mask.*t2meanimgs;
figure;
imagesc(t2_mean_masked); axis image;
set(gca,'XTick',[],'YTick',[])
set(gca,'fontsize',15)
j = jet;
j(1,:) = [0 0 0];
colorbar;
colormap(j);
ylabel(colorbar,'T_{2} mean (s)','fontsize',20);
r2_masked=sum(sum(r2))/nnz(r2)
mean_t2_masked=sum(sum(t2_mean_masked))/nnz(t2_mean_masked)

%%
ii=75; %110
jj=16; %41
modelfit = t2lamps(ii,jj).*exp(-t/t2limgs(ii,jj))+t2samps(ii,jj).*exp(-t/t2simgs(ii,jj));
decaycurve = squeeze(imgs(ii,jj,:));
figure;
plot(t,decaycurve,t,modelfit)

%%
% e = t(1:2);
% f = t(2:end);
% decaycurve1 = decaycurve(1:2);
% decaycurve2 = decaycurve(2:end);
% p1 = polyfit(e,log(decaycurve1'),1);
% p2 = polyfit(f,log(decaycurve2'),1);
% yfit1 = exp(polyval(p1,e));
% yfit2 = exp(polyval(p2,f));
% yfitoverall = zeros(20,1);
% yfitoverall(1:2,1) = yfit1;
% yfitoverall(3:20,1) = yfit2(2:19);
% figure;
% plot(t,decaycurve)
% set(gca, 'YScale', 'log');
% z = t';
% figure;
% semilogy(t,decaycurve,'*',t,(yfitoverall))
% figure;
% plot(t,yfitoverall)

%%
clims = [0 max(imgs(:))];
figure;imagescn(imgs,clims,[1 y(3)]); colormap(jet)

figure;
clims = [0 max(imgs(:))];
for i=1:y(3)
    subplot_tight(1,y(3)/1,i,0.002); imagesc(imgs(:,:,i),clims); colormap(gray);  axis image; set(gca,'XTick',[],'YTick',[]); colormap(jet);
end