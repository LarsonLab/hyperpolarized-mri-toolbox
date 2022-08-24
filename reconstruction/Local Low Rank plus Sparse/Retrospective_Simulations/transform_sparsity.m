close all; 
clc;
clear all;

%% Create Dataset
epi_kidney = 0;
if epi_kidney == 2
    load epi_rat_kidney.mat 

    epi_kidney = pa(:,:,[1 3 5 7 9 11 13 15 17 19 21 23 25]); %choose 13 timepoints
%     epi_kidney = pa(:,:,2:17);
    for i = 1:13
        epi_kidney_reformat(:,:,i) = imresize(epi_kidney(13:21,7:28,i)',[40 12]); %resize to 40x12x13
    end

    data = ifftnc_time(epi_kidney_reformat,3); %generate kspace
    data = data/max(data(:)); %normalize to 1
    data = imnoise(data,'gaussian',0,0.000008); %add noise to make kidney SNR ~50 and aorta SNR ~85
    pa = fftnc_time(data,3); %ground truth
    q = 0:5:60;
    qq = 0:4:60;
    for i = 1:size(pa,1)
        for j = 1:size(pa,2)
            hp001_spline(i,j,:) = spline(q,pa(i,j,:),qq);
        end
    end
    pa = hp001_spline;
elseif epi_kidney == 1
    load epi_TRAMP_tumor.mat

    for i = 1:16
        tumor1(:,:,i) = pa(:,:,i,8);
        tumor(:,:,i) = imresize(tumor1(3:12,3:16,i),[20,20]);
    end

    data = ifftnc_time(tumor,3);
    data = data/max(data(:));
%     data = imnoise(data,'gaussian',0,0.0000000005);
    pa = fftnc_time(data,3);
    
else
    load 'n15urea_3D_kspace_VFA.mat'

    load vasculature.mat
    vasculature = curve_original/curve_original(1,1);
    load left_kidney.mat
    kidney = curve_original/curve_original(1,1);

    images = fftnc_time(kspace_pf,4);

    original3D = images(40,:,:,:);
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
                    pa(i,j,k) = images(i,j,k)*kidney(1,k)*0.2;
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
    data = imnoise(data,'gaussian',0,0.00000005);
    pa = fftnc_time(data,3);
    
    q = 0:5:60;
    qq = 0:4:60;
    for i = 1:size(pa,1)
        for j = 1:size(pa,2)
            hp001_spline(i,j,:) = spline(q,pa(i,j,:),qq);
        end
    end
    pa = hp001_spline;
end
%%
y = size(pa);

pa_reformat = reshape(pa,[y(1)*y(2) y(3)]); %reshape into space-time matrix
pa = pa/max(pa(:)); %normalize to 1
num_coefficients = 10; %percent of coefficients to keep

%% Temporal FFT
transform = TempFFT(2)*pa_reformat; %transform to sparse domain (from k-t sparse-sense package on http://cai2r.net/resources/software/k-t-sparse-sense-matlab-code)
transform_3D_tempfft = reshape(transform,[y(1) y(2) y(3)]);
m = sort(abs(transform(:)),'descend'); %sort transform coefficients
ndx = floor(length(m)*num_coefficients/100); %set threholding level
thresh = m(ndx);
transform_thresh = transform.*(abs(transform) > thresh); %create thresholded transform
pa_sparse_tempfft = TempFFT(2)'*transform_thresh; %back to spatial domain
pa_sparse_tempfft = pa_sparse_tempfft/max(pa_sparse_tempfft(:));
recon_tempfft = reshape(pa_sparse_tempfft,[y(1) y(2) y(3)]); %back to 40x12x13 matrix
ssimval_tempfft = zeros(1,y(3)); ssimmap_tempfft = zeros(y(1),y(2),y(3));
for i = 1:y(3)
    [ssimval_tempfft(1,i),ssimmap_tempfft(:,:,i)] = ssim_index(abs(pa(:,:,i)),abs(recon_tempfft(:,:,i))); %structural similarity index
end
mean_ssimval_tempfft = mean(ssimval_tempfft);
diff_tempfft = pa-recon_tempfft; %difference image
rmse_tempfft = abs(RMSE(pa(:),recon_tempfft(:))/rms(pa(:))); %NRMSE
nRMSE_tempfft = norm(abs(recon_tempfft(:)) - abs(pa(:))) / norm(abs(pa(:)));

%% Wavelet along time
transform = Wavelet_1d('Daubechies',4,1,2)*pa_reformat; %transform to sparse domain (from sparseMRI package on http://people.eecs.berkeley.edu/~mlustig/Software.html)
transform_3D_wavelet = reshape(transform,[y(1) y(2) y(3)]);
m = sort(abs(transform(:)),'descend'); %sort transform coefficients
ndx = floor(length(m)*num_coefficients/100); %set threholding level
thresh = m(ndx);
transform_thresh = transform.*(abs(transform) > thresh); %create thresholded transform
pa_sparse_wavelet = Wavelet_1d('Daubechies',4,1,2)'*transform_thresh; %back to spatial domain
pa_sparse_wavelet = pa_sparse_wavelet/max(pa_sparse_wavelet(:));
recon_wavelet = reshape(pa_sparse_wavelet,[y(1) y(2) y(3)]); %back to 40x12x13 matrix
ssimval_wavelet = zeros(1,y(3)); ssimmap_wavelet = zeros(y(1),y(2),y(3));
for i = 1:y(3)
    [ssimval_wavelet(1,i),ssimmap_wavelet(:,:,i)] = ssim_index(abs(pa(:,:,i)),abs(recon_wavelet(:,:,i))); %structural similarity index
end
mean_ssimval_wavelet = mean(ssimval_wavelet);
diff_wavelet = pa-recon_wavelet; %difference image
rmse_wavelet = abs(RMSE(pa(:),recon_wavelet(:))/rms(pa(:))); %NRMSE
nRMSE_wavelet = norm(abs(recon_wavelet(:)) - abs(pa(:))) / norm(abs(pa(:)));

%% PCA
[coeff] = pca(pa_reformat);
transform = pa_reformat*coeff; %transform to sparse domain
transform_3D_pca = reshape(transform,[y(1) y(2) y(3)]);
m = sort(abs(transform(:)),'descend'); %sort transform coefficients
ndx = floor(length(m)*num_coefficients/100); %set threholding level
thresh = m(ndx);
transform_thresh = transform.*(abs(transform) > thresh); %create thresholded transform
pa_sparse_pca = transform_thresh*coeff'; %back to spatial domain
pa_sparse_pca = pa_sparse_pca/max(pa_sparse_pca(:));
recon_pca = reshape(pa_sparse_pca,[y(1) y(2) y(3)]); %back to 40x12x13 matrix
ssimval_pca = zeros(1,y(3)); ssimmap_pca = zeros(y(1),y(2),y(3));
for i = 1:y(3)
    [ssimval_pca(1,i),ssimmap_pca(:,:,i)] = ssim_index(abs(pa(:,:,i)),abs(recon_pca(:,:,i))); %structural similarity index
end
mean_ssimval_pca = mean(ssimval_pca);
diff_pca = pa-recon_pca; %difference image
rmse_pca = abs(RMSE(pa(:),recon_pca(:))/rms(pa(:))); %NRMSE
nRMSE_pca = norm(abs(recon_pca(:)) - abs(pa(:))) / norm(abs(pa(:)));

%% Total Variation along time
transform = TV_Temp_2nd_dim*pa_reformat; %transform to sparse domain
transform_3D_tv = reshape(transform,[y(1) y(2) y(3)]);
m = sort(abs(transform(:)),'descend'); %sort transform coefficients
ndx = floor(length(m)*num_coefficients/100); %set threholding level
thresh = m(ndx);
transform_thresh = transform.*(abs(transform) > thresh); %create thresholded transform
pa_sparse_tv = transform_thresh(:,1);
for i = 1:(y(3)-1)
    pa_sparse_tv(:,i+1) = pa_sparse_tv(:,i) + transform_thresh(:,i+1);
end
% pa_sparse_tv = TV_Temp_2nd_dim'*transform_thresh; %back to spatial domain
pa_sparse_tv = pa_sparse_tv/max(abs(pa_sparse_tv(:)));
recon_tv = reshape(pa_sparse_tv,[y(1) y(2) y(3)]); %back to 40x12x13 matrix
recon_tv = abs(recon_tv/max(abs(recon_tv(:))));
ssimval_tv = zeros(1,y(3)); ssimmap_tv = zeros(y(1),y(2),y(3));
for i = 1:y(3)
    [ssimval_tv(1,i),ssimmap_tv(:,:,i)] = ssim_index(abs(pa(:,:,i)),abs(recon_tv(:,:,i))); %structural similarity index
end
mean_ssimval_tv = mean(ssimval_tv);
diff_tv = pa-recon_tv; %difference image
rmse_tv = abs(RMSE(pa(:),recon_tv(:))/rms(pa(:))); %NRMSE
nRMSE_tv = norm(abs(recon_tv(:)) - abs(pa(:))) / norm(abs(pa(:)));
