function BiCsliceProfile
%This function will explore the impact of slice profile effects for a
%Gaussian RF pulse and three different flip schedules. You will also need
%the 'rf_tools' pulse design library for MATLAB to run this function.

%% Altering the flip angle schedule

%Gaussian pulse shape on the Varian scanner
load gauss.mat 
% rf = msinc(256,.5); %TBW = 2 sinc pulse

pH = 7.4;
pKa = 6.17;
AtoB = 10^(pH-pKa);
% kex = 0; %in 1/s
kex = 1.56; %in 1/s
% kex = 4.26; %in 1/s
% kex = 4.32; %in 1/s
% kex = 5.51; %in 1/s
TR = .067; %in s
% TR = .33; %in s
% TR = .9; %in s
nRF = 36;

% %Constant 2.78deg flip angle BiC, 25deg flip angle CO2, 64 RF pulses
% nRF = 64;
% flips(1,:) = 25/9*pi/180.*ones(1,nRF); %bic
% flips(2,:) = 25*pi/180.*ones(1,nRF); %co2

flips(1,:) = 25/9*pi/180.*ones(1,nRF); %bic
flips(2,:) = 25*pi/180.*ones(1,nRF); %co2

scale = ones(1,nRF);
% load gcorr.mat %to use gradient scaling saved prior
if nRF > length(scale) %fill in last values of scale with final value
    scale = [scale ones(1,nRF-length(scale))*scale(end)];
else
    scale = scale(1:nRF); %so that scale vector size equals flips size
end

% [cf, s, out] = slice_profiles_exch(rf,flips,AtoB,TR,kex);
[cf, s, out] = slice_profiles_exch(rf,flips,AtoB,TR,kex,scale);

% Calculate pH error (compare all excitations w/ 1st excitation)
%
pHerr = log10(sum(sum(abs(out.mxy_uncorr.a)))/sum(sum(abs(out.mxy_uncorr.b))) / ...
     (sum(abs(out.mxy_uncorr.a(:,1)))/sum(abs(out.mxy_uncorr.b(:,1)))));
disp(['Calculated pH error using all excitations: ' num2str(pHerr,3)]);

% % If you want to save the calculated gradient scaling, uncomment out below!
% %
% scale = mean([out.fwhm.a;out.fwhm.b]); %uses mean normalized FWHM for gradient correction
% save('gcorr.mat','scale');
% disp('Saved gradient scaling factor in gcorr.mat!')

% % Calculate pH at each number of excitations (cumulative)
% %
% pH = 

% %RF VFA flip schedule
% flips = [20.7 22.2 24.1 26.6 30 35.3 45 90] .* pi/180;
% [cf(2,:), s(2,:), gcf(2,:), out(2)] = slice_profiles(rf,flips);
% 
% %RF, ADC VFA flip schedule
% flips = [30.8 31.5 32.6 34.2 36.6 40.7 45 90] .* pi/180;
% [cf(3,:), s(3,:), gcf(3,:), out(3)] = slice_profiles(rf,flips);
%% Plotting

% Plot total BiC and CO2 signal (summed across slice) over all excitations
%
figure; 
leg=[{'SliceProf'},{'no-SliceProf'}];
subplot(121); plot(1:nRF,[out.signal_total.a' s(1,:)']); title('BiC signal'); 
legend(leg); axis([1 nRF 0 1]);
subplot(122); plot(1:nRF,[out.signal_total.b' s(2,:)']); title('CO2 signal'); legend(leg)
legend(leg); axis([1 nRF 0 1]);

%Simulate PSFs due to sampling (Note: if exchange is fast enough, BiC and 
% CO2 are broadened equally for a given pH)
[bpsf,bFWHM] = BiCpsf(out.signal_total.a',pH,flips(1,1),flips(2,1),1);
[cpsf,cFWHM] = BiCpsf(out.signal_total.b',pH,flips(1,1),flips(2,1),0);

% %set up legend
% for ii = 1:nRF
%     legendInfo{ii} = ['RF #' num2str(ii)];
% end
% clear ii
% 
% ctr = [1 4 7];
% h = figure(1234);
% for jj = 1:3
%     subplot(3,3,ctr(jj))
%     plot(abs(out(jj).mxy_uncorr))
%     xlabel('Position'), ylabel('Signal (a.u.)')
%     set(gca,'XTick',[])
%     legend(legendInfo)
%     ylim([0 0.6])
%     xlim([735 1265])
%     
%     subplot(3,3,ctr(jj)+1)
%     plot(abs(out(jj).mxy_gradcorr)) 
%     xlabel('Position'), ylabel('Signal (a.u.)')
%     set(gca,'XTick',[])
%     legend(legendInfo)
%     ylim([0 0.6])
%     xlim([735 1265])
%     
%     subplot(3,3,ctr(jj)+2)
%     %Simulated signal without gradient correction
%     plot(s(jj,:).*cf(jj,:),'color', [0.8500    0.3250    0.0980])
%     hold on
%     %Confirm that slice-select gradient corrected signal matches desired signal s
%     plot(sum(abs(out(jj).mxy_gradcorr),1)./max(max(sum(abs(out(jj).mxy_gradcorr),1))),'color', [0    0.4470    0.7410])
%     %Desired signal
%     plot(s(jj,:),'--r')    
%     xlim([1 8])
%     ylim([0 2.1])
%     xlabel('RF #')
%     ylabel('Normalized Signal')
%     legend('Simulated: No Correction','Simulated: Slice Correction','Desired','Location','south')
% end
% 
% display('Maximum signal deviation due to slice profile effects is:')
% display(['Constant 30deg: ' num2str(round(cf(1,end),2))])
% display(['RF VFA: ' num2str(round(cf(2,end),2))])
% display(['RF, ADC VFA: ' num2str(round(cf(3,end),2))])
% %% Altering the RF pulse shape
% %%It's important to note that this effect is a function of the RF pulse shape 
% %(excitation profile) and flip angle schedule. Altering either of these
% %will change the magnitude of the slice profile effect, as we'll explore
% %below
% 
% close all, clear all, clc
% load gauss.mat 
% nRF = 8;
% 
% %Constant 30deg flip angle, Gaussian RF
% flips = 30*pi/180.*ones(1,nRF);
% [cf, s, gcf, out] = slice_profiles(rf,flips);
% 
% %Constant 30deg flip angle, low TBW sinc
% rf = msinc(256,1);
% [cf(2,:), s(2,:), gcf(2,:), out(2)] = slice_profiles(rf,flips);
% 
% %Constant 30deg flip angle, high TBW sinc
% rf = msinc(256,4);
% g = 4.*ones(256,1);
% [cf(3,:), s(3,:), gcf(3,:), out(3)] = slice_profiles(rf,flips,g);
% 
% %% Plot the slice profiles and expected signal with and without gradient correction
% clc
% %set up legend
% for ii = 1:nRF
%     legendInfo{ii} = ['RF #' num2str(ii)];
% end
% clear ii
% 
% ctr = [1 4 7];
% h = figure(1234);
% for jj = 1:3
%     subplot(3,3,ctr(jj))
%     plot(abs(out(jj).mxy_uncorr))
%     xlabel('Position'), ylabel('Signal (a.u.)')
%     set(gca,'XTick',[])
%     legend(legendInfo)
%     ylim([0 0.6])
%     xlim([735 1265])
%     
%     subplot(3,3,ctr(jj)+1)
%     plot(abs(out(jj).mxy_gradcorr)) 
%     xlabel('Position'), ylabel('Signal (a.u.)')
%     set(gca,'XTick',[])
%     legend(legendInfo)
%     ylim([0 0.6])
%     xlim([735 1265])
%     
%     subplot(3,3,ctr(jj)+2)
%     plot(s(jj,:).*cf(jj,:),'color', [0.8500    0.3250    0.0980])
%     hold on
%     %Confirm that slice-select gradient corrected signal matches desired signal s
%     plot(sum(abs(out(jj).mxy_gradcorr),1)./max(max(sum(abs(out(jj).mxy_gradcorr),1))),'color', [0    0.4470    0.7410])
%     plot(s(jj,:),'--r')    
%     xlim([1 8])
%     ylim([0 1])
%     xlabel('RF #')
%     ylabel('Normalized Signal')
%     legend('Simulated: No Correction','Simulated: Slice Correction','Desired','Location','south')
% end
% 
% display('Maximum signal deviation due to slice profile effects is:')
% display(['Gaussian RF: ' num2str(round(cf(1,end),2))])
% display(['Low TBS Sinc: ' num2str(round(cf(2,end),2))])
% display(['High TBW Sinc: ' num2str(round(cf(3,end),2))])
end