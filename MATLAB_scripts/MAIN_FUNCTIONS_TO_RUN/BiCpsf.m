function [psf,FWHM] = BiCpsf(S,pH,FAb,FAc,bflag)
% BiCpsf.m - Simulates the point-spread function due to signal amplitude
% fluctuation during sampling of HP signal. Assumes centric encoding
% (subsequent excitations spiral out from center of k-space)
% 
% Updated 4/29/17 by D. Korenchan
%
% INPUTS:
%    S - matrix containing column vectors of signal amplitudes. Number of rows
%        must be the square of an integer
%    pH - vector containing pH value for each corresponding signal amplitude 
%         column vector in S
%    FAb - flip angle on bicarbonate (°)
%    FAc - flip angle on CO2 (°)
%    bflag - Tells script which metabolite it's simulating (1 is bicarb, 
%            2 is CO2) 
% OUTPUTS:
%    psf - 3D matrix containing all 2D PSF matrices from S (arrayed along
%    3rd dimension)
%    FWHM - vector containing FWHM of each corresponding PSF. Last value is
%    the FWHM of the PSF with equal k-space weighting

% Initialize variables
Npe = sqrt(size(S,1)); %square matrix dimension
Kwt = zeros(Npe,Npe,size(S,2)); %k-space weighting matrix
Nplot = 1000; %matrix dimension for plotting; must be odd/even if Npe is odd/even
Kwt_pad = zeros(Nplot,Nplot,size(S,2)+1); %zero-filled k-space weighting matrix
psf = Kwt_pad; %spatial PSF (FT of Kwt_pad)
if bflag %check which metabolite, adjust graph labels accordingly
    metname = 'BiC';
else
    metname = 'CO2';
end

% Normalize S
S = S ./ (ones(size(S,1),1) * max(S,[],1));

% Populate k-space sampling matrix, centric encoding
for m = 1:size(S,2)
    i = floor(Npe / 2 + 1); %initial row for populating k-weighting matrix
    j = floor(Npe / 2 + 1); %initial column for populating k-weighting matrix
    Kwt(i,j,m) = S(1,m);
    for k = 1:(Npe - 1)
        i = i + (-1) ^ k;
        Kwt(i,j,m) = S(1 + k^2,m);
        j = j + (-1) ^ k;
        j2 = j + (k - 1) * (-1)^k;
        Kwt(i,j:(-1) ^ (k):j2,m) = S(2 + k^2:2 + k^2 + (k - 1),m);
        j = j2;
        i = i - (-1) ^ k;
        i2 = i - (k-1) * (-1)^k;
        Kwt(i:(-1) ^ (k - 1):i2,j,m) = S((k + 1)^2 - (k - 1):(k + 1)^2,m);
        i = i2;
    end
end

% Prepare matrix for plotting, 2DFFT, plot
Npad = (Nplot - Npe) / 2;
for m = 1:size(S,2)
    Kwt_pad(Npad + 1:Npad + Npe,Npad + 1:Npad + Npe,m) = Kwt(:,:,m);
    psf(:,:,m) = 1 / (Npe ^ 2) * fftshift(fft2(fftshift(Kwt_pad(:,:,m))));
%    figure
%    surf(abs(psf(:,:,m)))
%    title(['PSF, pH = ' num2str(pH(m))]) 
end

% Plot PSF of equal weighting in all k-space points
Kwt_pad(Npad + 1:Npad + Npe,Npad + 1:Npad + Npe,end) = ones(Npe);
psf(:,:,end) = 1 / (Npe ^ 2) * fftshift(fft2(fftshift(Kwt_pad(:,:,end))));
figure
subplot(1,3,1)
surf(abs(psf(:,:,end)),'EdgeColor','none'); axis off;
title('PSF, no k-space weighting','FontSize',16)

% Determine FWHM of each PSF
for m = 1:size(psf,3)
    [psf_max,max_index] = max(max(abs(psf(:,:,m))));
    psf_search = abs(psf(:,max_index,m)) - psf_max / 2;
    [~,halfmax_index] = min(abs(psf_search));
    FWHM(m) = 2 * (max_index - halfmax_index);
end

% Plot PSFs with smallest and largest FWHMs and print ratio to uniform k-space PSF
[FWHM_max,FWHM_index] = max(FWHM(1:size(S,2)));
multmax = find(FWHM(1:size(S,2)) == FWHM_max);
[FWHM_min,FWHM_mindex] = min(FWHM(1:size(S,2)));
multmin = find(FWHM(1:size(S,2)) == FWHM_min);
subplot(1,3,2)
surf(abs(psf(:,:,FWHM_index)),'EdgeColor','none'); axis off
title([metname ' PSF, pH = ' num2str(pH(FWHM_index))],'FontSize',16) % ', BiC ' num2str(FAb) '{\circ} flip, CO_2 ' num2str(FAc) '{\circ} flip'
FWHM_ratio = FWHM_max / FWHM(end);
fprintf(1,['Max broadening: the ' metname ' PSF at pH ' num2str(pH(FWHM_index)) ' is ' num2str(FWHM_ratio,3) ' times broader than with uniform k-space weighting \n'])
fprintf(1,['(Max broadening found at ' num2str(length(multmax)) ' pH value(s)) \n'])
subplot(1,3,3)
surf(abs(psf(:,:,FWHM_mindex)),'EdgeColor','none'); axis off
title([metname ' PSF, pH = ' num2str(pH(FWHM_mindex))],'FontSize',16) % ', BiC ' num2str(FAb) '{\circ} flip, CO_2 ' num2str(FAc) '{\circ} flip'
FWHM_ratio = FWHM_min / FWHM(end);
fprintf(1,['Min broadening: the ' metname ' PSF at pH ' num2str(pH(FWHM_mindex)) ' is ' num2str(FWHM_ratio,3) ' times broader than with uniform k-space weighting \n'])
fprintf(1,['(Min broadening found at ' num2str(length(multmin)) ' pH value(s)) \n'])