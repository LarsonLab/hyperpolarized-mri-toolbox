%% iterative method: it is better to get lsqr work first!!!!
% gradient waveforms
G=[linspace(0,1,23),ones(1,454),linspace(1,0,23),linspace(0,-1,23),-ones(1,454),linspace(-1,0,23)];  
K=cumsum(G);
K=K/max(K)*2*pi-pi;

Kx=[K(38:10:470),K(538:10:970)];
Kx=Kx/max(Kx)*pi;

Kx1= [linspace(-pi,pi,50),linspace(pi,-pi,50)];
Kx = [Kx1(4:47),Kx1(54:97)];

Kxx=repmat(Kx,[1 25]);
Ky=linspace(-pi,pi,22);
[kx,ky]=meshgrid(Kxx,Ky);
kx=kx';
ky=ky';

Kf=linspace(-pi,pi,2500);
for ind=1:50
    Kf1((1+(ind-1)*44):(ind*44))= Kf((4+(ind-1)*50):(47+(ind-1)*50));
end
Kf=Kf1';

kf=repmat(Kf,[1 N*2]);

kxx=kx(:);
kyy=ky(:);
kff=kf(:);
%%
st_EPSI=nufft_init([kyy,kxx,kff], [N*2,N*2,50], [5,5,5], [4*N,4*N,100], [N,N,25],'minmax:kb');
%%
a=rawloadX('EPSI_1009_invivo');
a(:,1:2:end)=-a(:,1:2:end); % notice there is still chopping exsiting
%% off-isocenter shift correction ?optional?
% shift=7.3; % unit is mm from scanner
% G_unit=[linspace(0,1,3),ones(1,44),linspace(1,0,3)];
% G=[G_unit,-G_unit];
% gradient=repmat(G,[1 25]);
% Kx=cumsum(gradient);
% 
% Wk=max(Kx)-min(Kx);
% Kx=Kx-Wk/2;
% 
% Kx=Kx/max(Kx)*0.1389;
% 
% phase=repmat(Kx'*shift,[1 N*2]);
% 
% a=a.*exp(-1j*phase*2*pi);

%%
% get the sample from EPSI plautau
temp=zeros(44*50,22);
for ind=1:50
temp((1+(ind-1)*44):(ind*44),:)=a((4+(ind-1)*50):(47+(ind-1)*50),:);
end
data=temp;
data=data/length(data(:));

%% Interpolated and shifted version!
b=zeros(2500,22);
for ind=1:22
    x=1:10:25000;
    y=a(:,ind);
    x1=1:25000;
    y1=interp1(x,y,x1,'spline');
    b(:,ind)=y1(5:10:end);
end

temp=zeros(44*50,22);
for ind=1:50
temp((1+(ind-1)*44):(ind*44),:)=b((4+(ind-1)*50):(47+(ind-1)*50),:);
end
data=temp;
data=data/length(data(:));
%% iterative NUFFT version with regularization
lambda=0.5; % regularization term
%data=ones(size(data));
E = @(x,tr) NUFT(x,tr,st_EPSI);
F = @ (x,tr) NUFT_reg(x,tr,st_EPSI,lambda); % regularized iteration
spec_spatial = lsqr(E, data(:),1e-4,30);  % number of 30 is ok, 50 is noisy
spec_spatial = reshape(spec_spatial,st_EPSI.Nd);

%% The following part is to display

% display spatial-spectral info
figure,plot_voxels(abs(spec_spatial))

%% perform the shift
kdata=fftnc(spec_spatial);
% linear phase compensation in kspace
shift=7.3; % unit is mm

Kx=linspace(-1,1,22);
Kx=Kx/max(Kx)*0.1389;

phase=repmat(Kx*shift,[N*2 1]);

kdata1=zeros(22,22,50);
for ind=1:50
    kdata1(:,:,ind)=kdata(:,:,ind).*exp(-1j*2*pi*phase);
end

ss=ifftnc(kdata1);

%%
figure,plot_voxels(abs(ss))
ssp(:,:,:)=ss(:,:,end:-1:1);

%% extract 2D image
temp=zeros(44);

im_py=sum(abs(ssp(:,:,43:48)),3);


temp(12:33,12:33)=fft2c(im_py);%*exp(-1j*2*pi*phase.');

im_py=ifft2c(temp);

figure,imagesc(abs(im_py)),colormap(gray)

%% alanine
temp=zeros(44);

im_ala=sum(abs(ssp(:,:,23:28)),3);

temp(12:33,12:33)=fft2c(im_ala);

im_ala=ifft2c(temp);

%pyruvate hydrate
temp=zeros(44);

im_pyh=sum(abs(ssp(:,:,15:18)),3);

temp(12:33,12:33)=fft2c(im_pyh);

im_pyh=ifft2c(temp);

%lactate hydrate
temp=zeros(44);

im_lac=sum(abs(ssp(:,:,3:7)),3);

temp(12:33,12:33)=fft2c(im_lac);

im_lac=ifft2c(temp);

%%
figure,subplot(144),imagesc(abs(im_py)),colormap(gray)
subplot(143),imagesc(abs(im_ala)),colormap(gray)
subplot(142),imagesc(abs(im_pyh)),colormap(gray)
subplot(141),imagesc(abs(im_lac)),colormap(gray)
%% 2DFT image

k_2d=a(4:2:47,:);
k_2d=k_2d.'.*exp(-1j*2*pi*phase);
im_2d=ifft2c(k_2d);
%im_2d=im_2d.';

temp=zeros(44);

temp(12:33,12:33)=fft2c(im_2d);

im_2d=ifft2c(temp);
figure,imagesc(abs(im_2d)),colormap(gray)