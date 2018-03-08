close all;
clear all;

tic
for aa = 1:10
load 'n15urea_3D_kspace_VFA.mat'

load vasculature.mat
vasculature = curve_original/curve_original(1,1);
load left_kidney.mat
kidney = curve_original/curve_original(1,1);

images = fftnc_time(kspace_pf,4);
% images = image; 

original3D = images(40,:,:,:); %40 for paper
original3D = squeeze(original3D);
original3D = original3D/max(original3D(:));

noise = std2(original3D(36:40,1:3,1));
imgs = original3D/noise;

mask = (abs(imgs))>7;

images = repmat(original3D(:,:,3),[1 1 13]);

y = size(images);
r=rand(y(1),y(2),y(3));

for i=1:y(1)
    for j=1:y(2)
        for k=1:y(3)
            if mask(i,j,3)==1
                pa(i,j,k) = images(i,j,k)*kidney(1,k)*0.27;
            else
                pa(i,j,k) = images(i,j,k)+0.01*r(i,j,k);
            end
        end
    end
end

for i =1:y(3)
pa(18:20,6:8,i) = pa(18:20,6:8,i)*vasculature(1,i)*0.01;
pa(18:20,6:8,2) = pa(18:20,6:8,2)*1.3;
pa(18:20,6:8,3) = pa(18:20,6:8,3)*1.12;
end

noise2 = std2(pa(36:40,1:3,1));
imgs_pa = pa/noise2;

% for i=1:y(1)
%     for j=1:y(2)
%         for k=1:y(3)
%             if imgs_pa(i,j,k)>=7
%                 pa(i,j,k) = pa(i,j,k)*0.1;
%             else
%                 pa(i,j,k) = pa(i,j,k);
%             end
%         end
%     end
% end

data = fftnc_time(pa,3);
pa = fftnc_time(data,3);
data = data/max(data(:));
data = imnoise(data,'gaussian',0,0.00000001);
% save('retrospective_bSSFP_kidney_kspace','data');
y=size(data);
test = reshape(pa,[y(1)*y(2) y(3)]);
[u,v,t] = svd(test,0);
keep = 1:2;
glr_new = u(:,keep)*v(keep,keep)*t(:,keep)';
glr_new = reshape(glr_new,[y(1) y(2) y(3)]);
[u_llr,vllr,t_llr] = svd(reshape(pa(19:26,6:8,:),[24 13]));
[u_llr2,vllr2,t_llr2] = svd(reshape(pa(27:34,6:8,:),[24 13]));
v = v/max(v(:));
vllr = vllr/max(vllr(:));
vllr2 = vllr2/max(vllr2(:));
x=1:13;
figure;semilogy(x,max(v,[],2),x,max(vllr2(1:13,:),[],2),'*-'); ylim([10^(-5) 1]); xlim([1 9]);

%% Undersampling Pattern
Sampling_Percentage = 25;
str = ['pattern_' int2str(Sampling_Percentage) '_percent_rat.mat'];
load (str);

kdata = data.*pattern(:,:,:);

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
mask = abs(kdata)>0;
    
DN = [y(1) y(2) y(3)];
timedim = 3; 
wavelet_time_scale = 1;
sparse = 0;
%% L+S reconstruction **************

if sparse == 0;
    
    for i=1:1 %10 used for paper
    param.E=pnDFT_time(mask,DN,timedim,1,2);
    param.d=kdata(:,:,:);
    % param.T=Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);
    param.T= TempFFT(3);
    param.lambda_L=0.01; %0.01
    param.lambda_S=0.001; %0.001
    param.nite=2500; 
    param.tol=1.5e-3; %2e-3 for llrps; 4e-3 for LLR
    param.bsize = [8,3];
    param.mu = 2;

    [L,S] = lps_ist_2D_LLR_plus_S_PCA_V2(param); 

    LplusS(:,:,:,i)=L+S;
   
    end
else
    param.E=pnDFT_time(undersample,DN,timedim,1,2);
    param.W = TempFFT(3);param.L1Weight=0.01;
    param.TV = TVOP_3D();param.TVWeight=0.001;
    param.y = kdata;
    param.nite = 48;
    param.display=1;

    % linear reconstruction
    recon_dft=param.E'*kdata;

    % k-t SPARSE-SENSE reconstruction
    fprintf('\n k-t SPARSE-SENSE reconstruction \n')
    
    recon_cs=recon_dft;
    % the non-linear CG algorithm will re-start two times to improve
    % convergence using the last result as starting point
    for h=1:8,
%         recon_cs = CSL1NlCg(recon_cs,param);
        recon_cs = CSL1NlCg_pca(recon_cs,param);
    end
    LplusS = recon_cs;
    
end

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
ssim_map = zeros(y(1),y(2),y(3));
for i=1:y(3)
    [ssimval(1,i), ssim_map(:,:,i)] = ssim_index(abs(LplusS(:,:,i)),abs(pa(:,:,i)));
end
mean_ssimval = mean(ssimval(1,1:end));
nRMSE = norm(LplusS(:) - pa(:)) / norm(pa(:));
nRMSE_abs = norm(abs(LplusS(:)) - abs(pa(:))) / norm(abs(pa(:)));
nRMSE_abs_dft = norm(abs(recon_dft(:)) - abs(pa(:))) / norm(abs(pa(:)));
%%
vasculature_LplusS = zeros(1,y(3));
vasculature_orginal = zeros(1,y(3));
for i=1:y(3)
    vasculature_LplusS(1,i) = mean2(abs(LplusS(21:23,6:8,i)));
    vasculature_original(1,i) = mean2(abs(pa(21:23,6:8,i)));
    vasculature_dft(1,i) = mean2(abs(recon_dft(21:23,6:8,i)));
end
x = 1:y(3);
figure;
plot(x,vasculature_original,x,vasculature_LplusS,'--',x,vasculature_dft,'LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('Ground Truth','LLR+S','Zero-Filled Recon','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')
%%
kidney_LplusS = zeros(1,y(3));
kidney_orginal = zeros(1,y(3));
for i=1:y(3)
    kidney_LplusS(1,i) = mean2(abs(LplusS(26:36,1:7,i)));
    kidney_original(1,i) = mean2(abs(pa(26:36,1:7,i)));
    kidney_dft(1,i) = mean2(abs(recon_dft(26:36,1:7,i)));
end
x = 1:y(3);
figure;
plot(x,kidney_original,x,kidney_LplusS,'--',x,kidney_dft,'LineWidth',2)
xlabel('Timepoints','FontSize',20,'FontWeight','bold')
ylabel('Signal','FontSize',20,'FontWeight','bold') 
legend('Ground Truth','LLR+S','Zero-Filled Recon','Location','northeast')
set(gca,'FontSize',16,'FontWeight','bold')

%%
figure;
for i=1:y(3)
    subplot_tight(1,y(3),i,0.001); imagesc(rot90(abs(LplusS(:,:,i))),[0 1]); colormap(jet);  axis image; set(gca, 'XTick', [], 'YTick', [])
end

figure;
for i=1:y(3)
    subplot_tight(1,y(3),i,0.001); imagesc(rot90(abs(pa(:,:,i))),[0 1]); colormap(jet);  axis image; set(gca, 'XTick', [], 'YTick', [])
end

figure;
for i=1:y(3)
    subplot_tight(1,y(3),i,0.001); imagesc(rot90(abs(differ(:,:,i))),[0 1]); colormap(jet);  axis image; set(gca, 'XTick', [], 'YTick', [])
end

figure;
for i=1:y(3)
    subplot_tight(1,y(3),i,0.001); imagesc(rot90(abs(ssim_map(:,:,i))),[0 1]); colormap(jet);  axis image; set(gca, 'XTick', [], 'YTick', [])
end

if aa <10
	close all;
end
ssimval2(1,aa) = mean_ssimval;
nRMSE2(1,aa) = nrmse2;
aa
end
mean_ssimval2 = mean(ssimval2)
mean_nRMSE2 = mean(nRMSE2)
toc