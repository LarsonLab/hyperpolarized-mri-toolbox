%Example recon of 3D dynamic bSSFP data with LLR+S recon using LLR and PCA. Using a
%parfor loop to speed up recosntruction
%Rat abdomen (kidney) prospective data 
%Eugene Milshteyn 08/02/2016

load kspace_c2pyr_prospective %C2-pyruvate prospective rat data
% load kspace_lactate_prospective %C1-lactate prospective rat data
% load kspace_n5urea_prospective %N15Urea prospective rat data

close all; clear all;
tic
y = size(kspace);
dyn = y(4);
data = fft(kspace,[],1);
data = circshift(data,[y(1)/2 0 0 0]);
data = data/max(data(:));
% data = kspace;
undersample = zeros(y(1),y(2),y(3),dyn);
for i=1:y(1)
    for j=1:y(2)
        for k=1:y(3)
            for l=1:dyn
                if abs(kspace(i,j,k,l))>0
                    undersample(i,j,k,l) = 1;
                elseif abs(kspace(i,j,k,l)) == 0
                    undersample(i,j,k,l)=0;
                end
            end
        end
    end
end

data = permute(data,[2 3 4 1]);
undersample = permute(undersample,[2 3 4 1]);
DN = [y(2) y(3) dyn];
timedim = 3; 
wavelet_time_scale = 1;
LplusS = zeros(y(2),y(3),y(4),y(1),iter);
L = zeros(y(2),y(3),y(4),y(1),iter);
S = zeros(y(2),y(3),y(4),y(1),iter);
sparse = 0; %0 Means do LLR+S recon; 1 means do just sparse recon (nonlin conjugate gradient from NYU Group)
lambda_L=0.01;
lambda_S=0.001;
nite=2500;
tol=1.5e-3;
bsize = [8,3];
mu = 2;
%% L+S reconstruction ******************************************************
if sparse == 0  
    parfor i = 1:y(1)
        currentpattern = squeeze(undersample(:,:,:,i));
        E=pnDFT_time(currentpattern,DN,timedim,1,2);
        currentslice = squeeze(data(:,:,:,i));     
        [L,S] = LLR_plus_S_PCA_parfor(currentslice,E,lambda_L,lambda_S,nite,tol,bsize,mu);
        L(:,:,:,i) = L; %Local low rank portion
        S(:,:,:,i) = S; %Sparse portion
        LplusS(:,:,:,i)=L(:,:,:,i)+S(:,:,:,i); %LLR+S
    end 
else
    for i=1:y(1)
    param.E=pnDFT_time(undersample,DN,timedim,1,2);
    param.W = TempFFT(3);
    %     param.W=Wavelet_1d('Daubechies',4,wavelet_time_scale,timedim);
    param.L1Weight=0.05;
    param.TV = TVOP_3D();param.TVWeight=0.01;
    param.y = kdata(:,:,:,i);
    param.nite = 48;
    param.display=1;

    % linear reconstruction
    recon_dft(:,:,:,i)=param.E'*kdata(:,:,:,i);

    % k-t SPARSE-SENSE reconstruction
    fprintf('\n k-t SPARSE-SENSE reconstruction \n')
    tic
    recon_cs(:,:,:,i)=recon_dft(:,:,:,i);
    % the non-linear CG algorithm will re-start two times to improve
    % convergence using the last result as starting point
    for h=1:3,
        recon_cs(:,:,:,i) = CSL1NlCg(recon_cs(:,:,:,i),param);
    end
    LplusS(:,:,:,i)=recon_cs(:,:,:,i); %Sparse recon only; kept LplusS name for rest of script
    end
end

%% Correct for Voxel Shifts
for i = 1:size(LplusS,3)
    for j = 1:size(LplusS,4)
        LplusS(:,:,i,j) = flipud(LplusS(:,:,i,j));
        LplusS(:,:,i,j) = circshift(LplusS(:,:,i,j),[1 -1]);
        L(:,:,i,j) = flipud(L(:,:,i,j));
        L(:,:,i,j) = circshift(L(:,:,i,j),[1 -1]);
        S(:,:,i,j) = flipud(S(:,:,i,j));
        S(:,:,i,j) = circshift(S(:,:,i,j),[1 -1]);
    end
end


%Permute into [freq_encode phase_encode1 phase_encode2 timepoints]
LplusS = permute(LplusS,[4 1 2 3]); 
L = permute(L,[4 1 2 3]);
S = permute(S,[4 1 2 3]);

LplusS2 = zeros(y(1),y(2),y(3),y(4));
L2 = zeros(y(1),y(2),y(3),y(4));
S2 = zeros(y(1),y(2),y(3),y(4));

s = y(3):-1:1;

%Reorder Data to match intended acquisition: reorder slices and flip 2nd
%dimension
for i = 1:y(3)
    for j = 1:y(4)
        LplusS2(:,:,i,j) = fliplr(LplusS(:,:,s(i),j));
        L2(:,:,i,j) = fliplr(L(:,:,s(i),j));
        S2(:,:,i,j) = fliplr(S(:,:,s(i),j));
    end
end

LplusS = LplusS2;
L = L2;
S = S2;

toc

%% Save files
% save('LplusS_C1pyr_20170530_exp1_iter1_lpsPCA_V1','LplusS');
% save('L_C1pyr_20170530_exp1_iter1_lpsPCA_V1','L');
% save('S_C1pyr_20170530_exp1_iter1_lpsPCA_V1','S');

%% SNR 
imgnoise = std2(LplusS(1:10,1:10,1,end));
LplusS_SNR = LplusS/imgnoise;
