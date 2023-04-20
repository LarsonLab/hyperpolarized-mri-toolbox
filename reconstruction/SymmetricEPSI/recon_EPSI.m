function [im_zp Spect_data] = recon_EPSI(input);
% Nyquist EPSI: written by Wenwen on 05022014
% 250 *100 (by 10)
% each unit: 54 ramp up, 54 ramp down, 142 on plat
x=rawloadX(input);
x = x(15:2514,:);
x(:,1:2:end) = -x(:,1:2:end); % chopping effect
[N_rd N] = size(x);
Gx=([linspace(0,1,54),ones(1,142),linspace(1,0,54)]);
Kx = cumsum(Gx);
Kx = Kx/max(Kx) - 0.5;
Kx = Kx(1:10:end);
Ky = linspace(-0.5,0.5,N);
[kx ky] =meshgrid(Kx,Ky);
kx = kx.';
ky = ky.';

kxx = kx(:)*pi*2;
kyy = ky(:)*pi*2;
% approximate density compensation with gradient waveform
dx = repmat(Gx(1:10:end)',[1 N]);
dx = dx/max(dx(:));
% just ti display the image at the first frame
st_2D=nufft_init([-kyy,kxx], [N,N], [5,5], [2*N,2*N], [N/2,N/2],'minmax:kb');
% the beginning 
tmp=x(1:25,:).*dx;

image_2d = nufft_adj(tmp(:), st_2D);
ksp_zp = zeros(N*2);
ksp_zp(12:33,12:33)= fft2c(image_2d);
im_zp = ifft2c(ksp_zp);
figure,imagesc(abs(im_zp))
colormap(gray)


Kx = [Kx,Kx(end:-1:1)];
kkx = repmat(Kx,[1 50]);
kky = Ky;
kf=linspace(-pi,pi,N_rd);  %% should be -pi to pi, correct!!!
kf=kf';
kf = repmat(kf,[1 N]); 



[kx ky] =meshgrid(kkx,kky);
kx = kx.'*pi*2;
ky = ky.'*pi*2;


kxx=kx(:);
kyy=ky(:);
kff=kf(:);

% new version 090613: read the waveforms from the scanner!!! Standard!%%%
% as long as the Kx,Ky is correct

% this is for proton!!!
st=nufft_init([-kyy,kxx,kff], [N,N,50], [5,5,5], [2*N,2*N,50*2], [N/2,N/2,25],'minmax:kb');  % build the structure, remember to shift a little bit
%should be kyy,-kxx 
weight=repmat(dx,[100,1]);
weight=weight/sum(weight(:));
k_data=x.*weight;
k_data_new=k_data(:);  % shoule be the same size!!!!
Spect_data = nufft_adj(k_data_new, st);

im = sum(abs(Spect_data),3);
figure,imagesc(abs(im))
colormap(gray)

