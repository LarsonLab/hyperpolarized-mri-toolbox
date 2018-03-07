close all; 
clear all;

%% Create Dataset
epi_kidney = 3;
if epi_kidney == 3
    load t2mapping_lactate_retrospective_data
    pa = fftnc_time(data,3);
    bsize = [14,7];
elseif epi_kidney == 2
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
    bsize = [8,3];
elseif epi_kidney == 1
    load epi_TRAMP_tumor.mat

    for i = 1:16
        tumor1(:,:,i) = pa(:,:,i,8);
        tumor(:,:,i) = imresize(tumor1(3:12,3:16,i),[20,20]);
    end

    data = ifftnc_time(tumor,3);
    data = data/max(data(:));
    % data = imnoise(data,'gaussian',0,0.0000000005);
    pa = fftnc_time(data,3);
    bsize = [5,5];
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
    bsize = [8,3];
end

% load c13_kidney_eugene.mat %Jeremy EPI acquisition CFA
% 
% % epi_kidney = pa(:,:,[1 3 5 7 9 11 13 15 17 19 21 23 25]); %choose 13 timepoints
% epi_kidney = pa(:,:,2:17);
% for i = 1:16
%     epi_kidney_reformat(:,:,i) = imresize(epi_kidney(13:21,7:28,i)',[40 12]); %resize to 40x12x13
% end
% 
% data = ifftnc_time(epi_kidney_reformat,3); %generate kspace
% data = data/max(data(:)); %normalize to 1
% data = imnoise(data,'gaussian',0,0.000008); %add noise to make kidney SNR ~50 and aorta SNR ~85
% pa = fftnc_time(data,3); %ground truth
pa = pa/max(pa(:));
y = size(pa);
for ii = 2:y(3)
C = ii; 
compression_global = (y(1)*y(2)*y(3))/(C*(y(1)*y(2)+y(3)-C)) %compression ratio
compression_ratio(1,ii-1) = compression_global;
%% GLR
glr = reshape(pa,[y(1)*y(2) y(3)]);
[U,V,T] = svd(glr,0);
keep = 1:floor(C);
glr_new = U(:,keep)*V(keep,keep)*T(:,keep)';
glr_new = reshape(glr_new,[y(1),y(2),y(3)]);

pa = pa/max(pa(:));
glr_new = glr_new/max(glr_new(:));

%% Sparse Only
% sparse_transform = TempFFT(2);
%     sparse_transform = Wavelet_1d('Daubechies',4,1,2);
[coeff] = pca(glr);
transform = glr*coeff;

m = sort(abs(transform(:)),'descend');
ndx = floor((y(1)*y(2)*y(3))/compression_global);
thresh = m(ndx);
transform_thresh = transform.*(abs(transform) > thresh);
slice_denoise = transform_thresh*coeff';

slice_denoise = slice_denoise/max(slice_denoise(:));
slice_denoise = reshape(slice_denoise,[y(1),y(2),y(3)]);

%% GLR plus Sparse
C_rank = C-1;
glr = reshape(pa,[y(1)*y(2) y(3)]);
[U,V,T] = svd(glr,0);
keep = 1:floor(C_rank);
glr_new2 = U(:,keep)*V(keep,keep)*T(:,keep)';
glr_new2 = reshape(glr_new2,[y(1),y(2),y(3)]);

% glr_new = glr_new/max(glr_new(:));
compression = compression_global;
n = y(1)*y(2)*y(3);
nL = C_rank*(y(1)*y(2)+y(3)-C_rank);
sparse_matrix = pa-glr_new2;


slice = reshape(sparse_matrix,[y(1)*y(2) y(3)]);
slice = squeeze(slice);


% sparse_transform = TempFFT(3);
%     sparse_transform = Wavelet_1d('Daubechies',8,1,2);
sparse_matrix = reshape(sparse_matrix,[y(1)*y(2) y(3)]);
coeff = pca(sparse_matrix);
transform = sparse_matrix*coeff;
m = sort(abs(transform(:)),'descend');
% ndx = floor(((y(1)*y(2)*y(3))/compression)-(C_rank*(y(1)*y(2)+y(3)-C_rank)));
ndx = floor(n/compression-nL);
thresh = m(ndx);
transform_thresh = transform.*(abs(transform) > thresh);
slice_LplusS = transform_thresh*coeff';
% slice_LplusS = slice_LplusS/max(slice_LplusS(:));
slice_LplusS = reshape(slice_LplusS,[y(1),y(2),y(3)]);


LplusS = glr_new2+slice_LplusS;
LplusS = LplusS/max(LplusS(:));
% pa = pa/max(pa(:));

%% LLR plus Sparse
% if epi_kidney == 0||2
%     bsize = [8,3];
% elseif epi_kidney == 1
%     bsize = [5,5];
% else
%     bsize = [14,7];
% end
zpadx = y(1)/bsize(1);
zpady = y(2)/bsize(2);
Y = pa;
nt = y(3);
randX = randperm(zpadx);
randY = randperm(zpady);
% Y = circshift(Y, [randX(1)-1, randY(1)-1,  0]);
for blockx = 1:zpadx
    for blocky = 1:zpady

        tempb = Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:);
        tempb = reshape(tempb,[bsize(1)*bsize(2),nt]);
        [Ut,Vt,St] = svd(tempb,0);
        keep = 1:floor(C_rank);
        tempb = Ut(:,keep)*(Vt(keep,keep))*St(:,keep)';
        tempb = reshape(tempb,[bsize(1),bsize(2),nt]);
        Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:) = tempb;

    end
end
% Y = circshift(Y, [-randX(1)+1, -randY(1)+1,  0]);
glr_new_LLR = Y;


sparse_matrix2 = pa-glr_new_LLR;


slice2 = reshape(sparse_matrix2,[y(1)*y(2) y(3)]);
slice2 = squeeze(slice2);

% sparse_transform = TempFFT(3);
%     sparse_transform = Wavelet_1d('Daubechies',8,1,2);
sparse_matrix2 = reshape(sparse_matrix2, [y(1)*y(2) y(3)]);
coeff = pca(sparse_matrix2);
transform2 = sparse_matrix2*coeff;
m = sort(abs(transform2(:)),'descend');
% ndx = floor(((y(1)*y(2)*y(3))/compression)-(C_rank*(y(1)*y(2)+y(3)-C_rank)));
ndx = floor(n/compression-nL);
thresh = m(ndx);
transform_thresh2 = transform2.*(abs(transform2) > thresh);
slice_LplusS2 = transform_thresh2*coeff';
% slice_LplusS = slice_LplusS/max(slice_LplusS(:));
slice_LplusS2 = reshape(slice_LplusS2,[y(1),y(2),y(3)]);


LplusS2 = glr_new_LLR+slice_LplusS2;
LplusS2 = LplusS2/max(LplusS2(:));
% pa = pa/max(pa(:));

%% LLR
% if epi_kidney == 0||2
%     bsize = [8,3];
% elseif epi_kidney == 1
%     bsize = [5,5];
% else
%     bsize = [14,7];
% end
zpadx = y(1)/bsize(1);
zpady = y(2)/bsize(2);
Y = pa;
nt = y(3);
randX = randperm(zpadx);
randY = randperm(zpady);
% Y = circshift(Y, [randX(1)-1, randY(1)-1,  0]);
for blockx = 1:zpadx
    for blocky = 1:zpady

        tempb = Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:);
        tempb = reshape(tempb,[bsize(1)*bsize(2),nt]);
        [Ut,Vt,St] = svd(tempb,0);
        keep = 1:floor(C);
        tempb = Ut(:,keep)*(Vt(keep,keep))*St(:,keep)';
        tempb = reshape(tempb,[bsize(1),bsize(2),nt]);
        Y((blockx-1)*bsize(1)+1:bsize(1)*blockx,(blocky-1)*bsize(2)+1:bsize(2)*blocky,:) = tempb;

    end
end
% Y = circshift(Y, [-randX(1)+1, -randY(1)+1,  0]);
glr_new_LLR = Y;

%% Calculate rmse and differences
diff_glr = pa-glr_new;
rmse_glr = RMSE(pa(:),glr_new(:));
for i = 1:y(3)
    [mssim_glr(1,i),ssim_map_glr(:,:,i)] = ssim_index(pa(:,:,i),glr_new(:,:,i));
    rmse_test_glr(:,i) = RMSE(pa(:,:,i),glr_new(:,:,i));
end
mean_mssim_glr = abs(mean(mssim_glr(1,:)));
mean_rmse_test_glr = abs(mean(rmse_test_glr(:))/rms(pa(:)));
nRMSE_glr(1,ii) = norm(abs(glr_new(:)) - abs(pa(:))) / norm(abs(pa(:)));


diff_slice_denoise = pa-slice_denoise;
rmse_sparse = RMSE(pa(:),slice_denoise(:));
for i = 1:y(3)
    [mssim_sparse(1,i),ssim_map_sparse(:,:,i)] = ssim_index(pa(:,:,i),slice_denoise(:,:,i));
    rmse_test_sparse(:,i) = RMSE(pa(:,:,i),slice_denoise(:,:,i));
end
mean_mssim_sparse = abs(mean(mssim_sparse(1,:)));
mean_rmse_test_sparse = abs(mean(rmse_test_sparse(:))/rms(pa(:)));
nRMSE_sparse(1,ii) = norm(abs(slice_denoise(:)) - abs(pa(:))) / norm(abs(pa(:)));


diff_LplusS = pa-LplusS;
rmse_LplusS = RMSE(pa(:),LplusS(:));
for i = 1:y(3)
    [mssim_LplusS(1,i),ssim_map_LplusS(:,:,i)] = ssim_index(pa(:,:,i),LplusS(:,:,i));
    rmse_test_LplusS(:,i) = RMSE(pa(:,:,i),LplusS(:,:,i));
end
mean_mssim_LplusS = abs(mean(mssim_LplusS(1,:)));
mean_rmse_test_LplusS = abs(mean(rmse_test_LplusS(:))/rms(pa(:)));
nRMSE_LplusS(1,ii) = norm(abs(LplusS(:)) - abs(pa(:))) / norm(abs(pa(:)));


diff_LLR = pa-glr_new_LLR;
rmse_LLR = RMSE(pa(:),glr_new_LLR(:));
for i = 1:y(3)
    [mssim_LLR(1,i),ssim_map_LLR(:,:,i)] = ssim_index(pa(:,:,i),glr_new_LLR(:,:,i));
    rmse_test_LLR(:,i) = RMSE(pa(:,:,i),glr_new_LLR(:,:,i));
end
mean_mssim_LLR = abs(mean(mssim_LLR(1,3:end)));
mean_rmse_test_LLR = abs(mean(rmse_test_LLR(:))/rms(pa(:)));
nRMSE_LLR(1,ii) = norm(abs(glr_new_LLR(:)) - abs(pa(:))) / norm(abs(pa(:)));


diff_LplusS2 = pa-LplusS2;
rmse_LplusS2 = RMSE(pa(:),LplusS2(:));
for i = 1:y(3)
    [mssim_LplusS2(1,i),ssim_map_LplusS2(:,:,i)] = ssim_index(pa(:,:,i),LplusS2(:,:,i));
    rmse_test_LplusS2(:,i) = RMSE(pa(:,:,i),LplusS2(:,:,i));
end
mean_mssim_LplusS2 = abs(mean(mssim_LplusS2(1,3:end)));
mean_rmse_test_LplusS2 = abs(mean(rmse_test_LplusS2(:))/rms(pa(:)));
nRMSE_LplusS2(1,ii) = norm(abs(LplusS2(:)) - abs(pa(:))) / norm(abs(pa(:)));

end
%%
% 9:25,6:10
nRMSE_glr = fliplr(nRMSE_glr);
nRMSE_sparse = fliplr(nRMSE_sparse);
nRMSE_LplusS = fliplr(nRMSE_LplusS);
nRMSE_LLR = fliplr(nRMSE_LLR);
nRMSE_LplusS2 = fliplr(nRMSE_LplusS2);

for i = 1:y(3)
    dynamic_pa(1,i) = mean2(abs(pa(:,:,i)));
    dynamic_glr(1,i) = mean2(abs(glr_new(:,:,i)));
    dynamic_sparse(1,i) = mean2(abs(slice_denoise(:,:,i)));
    dynamic_LplusS(1,i) = mean2(abs(LplusS(:,:,i)));
    dynamic_LLR(1,i) = mean2(abs(glr_new_LLR(:,:,i)));
    dynamic_LplusS2(1,i) = mean2(abs(LplusS2(:,:,i)));
end
% x = 1:y(3)-1;
% figure;
% plot(x,dynamic_pa,'b',x,dynamic_glr,'r',x,dynamic_sparse,'g',x,dynamic_LplusS,'y',x,dynamic_LLR,'c',x,dynamic_LplusS2,'m')
%     
% figure;
% plot(x,dynamic_pa,'b',x,dynamic_glr,'r',x,dynamic_LplusS2,'m')
x = fliplr(compression_ratio);

figure;
plot(x,nRMSE_glr(1,1:y(3)-1),'b',x,nRMSE_sparse(1,1:y(3)-1),'g',x,nRMSE_LplusS2(1,1:y(3)-1),'k',x,nRMSE_LLR(1,1:y(3)-1),'r',x,nRMSE_LplusS(1,1:y(3)-1),'y','LineWidth',2)
xlabel('Compression Ratio','FontSize',20,'FontWeight','bold')
ylabel('nRMSE','FontSize',20,'FontWeight','bold') 
legend('Global Low Rank','Sparse Only','Local Low Rank plus Sparse','Local Low Rank','Low Rank plus Sparse','Location','NorthWest')
set(gca,'FontSize',16,'FontWeight','bold')
