function BiCexchEnhance
%BiCsim3.m
%
%Takes two flip angles for bicarb and CO2 and computes the z-magnetization 
%and signal magnitudes over a range of total excitations. Pulses are 
%applied to each resonance simultaneously. A Bloch-McConnell simulation is
%used inbetween excitations. The pH is calculated by summing 
%the signals, correcting for flip angle, and using the signal ratio in the
%Henderson-Hasselbalch equation. The script then plots the signal for each
%excitation, the total pH error as a function of true pH, and the CO2 
%signal enhancement over the no-exchange case as a function of true pH. The 
%script finally calls BiCpsf.m to plot/calculate each PSF and FWHM.
%
%Assumptions:
%   1. Bicarb and CO2 have equal T1 values due to rapid exchange
%
%INPUTS:
%   FAc, FAb:   Flip angle for each pulse on CO2, BiC (degrees)
%   Nexc:       Number of excitations (unitless)
%   pH:         Vector containing pH values of imaging volume (unitless)
%
%OUTPUTS:   
%   Sc, Sb:     Contains CO2, BiC signal magnitude produced by each flip 
%               (arb. units)
%   pH_cal:     Contains pH values calculated from BiC/CO2 signal
%               magnitudes over simulated pH range and #'s of summed
%               excitations (unitless)
%   ex_enhfac:  Enhancement factor of CO2 signal over the case with no 
%               exchange (ie. kex = 0) at each true pH value
%   psf:        Contains calculated point-spread functions for each pH
%               value and for a uniform k-space weighting
%   FWHM:       Vector of full-widths half-max for each PSF in psf
%
% 3/25/19:   Changed to specify timestep between points, rather than # of 
% points, for Bloch-McConnell simulation

%Initialize variables
%
FAc = 25; %CO2 tip angle (°)
FAb = FAc / 9; %BiC tip angle (°)
Nexc = 64; %total # of excitations
pH = 6.4:0.1:7.6; %vector of pH values to simulate
% kex = 0;
% kex = .168; %exchange rate between BiC and CO2 (1/s) - in free solution, no CA
% kex = 1.3; %exchange rate between BiC and CO2 (1/s) - measured in TRAMP tumors
kex = 5.51; %exchange rate between BiC and CO2 (1/s) - with CA at 7.55 ug/mL

TR = .067; %repetition time between excitations (s)
T1 = 3000000000000000000; %T1 of BiC and CO2 (s) (assumed equal)
pKa = 6.17; %pKa of BiC-CO2 in vivo at 37 degC
N = 1:Nexc;

dt = 0.01; %timestep between points for Bloch-McConnell (s)
np = ceil(TR / dt) + 1; %# of points to simulate over TR for exchange

%Initialize outermost loop, which will iterate each pH value
%
Mb = zeros(Nexc,length(pH));
Mc = Mb;
Mbtoc = Mc; %magnetization that transferred over from BiC to CO2 over TR
Sbtoc = Mbtoc; %signal that arose at end of TR from BiC->CO2 exchanged magnetization
Sb = Mb;
Sc = Sb;
Sbsum = Sb;
Scsum = Sbsum;
Mbtocsum = Scsum; %magnetization that transferred over from BiC to CO2 over each # of TRs
Sbtocsum = Mbtocsum; %total signal arising from BiC -> CO2 exchanged magnetization
pH_cal = Sb;
for m=1:size(pH,2)
    %Calculate ratio of bicarb:CO2
    BtoC=10 ^ (pH(m) - pKa);
    %Initial values of Mb and Mc are calculated using the
    %specified pH and the Henderson-Hasselbalch equation. Total equals 1
    Mb(1,m) = 1 / (1 / BtoC + 1);
    Mc(1,m) = 1 / BtoC / (1 / BtoC + 1);
    Mbtoc(1,m) = 0;
    Sbtoc(1,m) = 0;
    %Calculate kbc and kcb
    kbc = kex * Mc(1,m);
    kcb = kex * Mb(1,m);
    %Initialize loop, which will calculate the z-magnetization and the signal
    %magnitudes obtained after each pulse
    for i = 1:N(end)
        %Obtain signal (xy) magnetization for bicarb and CO2 on ith pulse
        Sc(i,m) = Mc(i,m) * sind(FAc);
        Sbtoc(i,m) = Mbtoc(i,m) * sind(FAc);
        Sb(i,m) = Mb(i,m) * sind(FAb);
        %Adjust bicarb and CO2 z-magnetizations after ith pulse
        Mc(i+1,m) = Mc(i,m) * cosd(FAc);
        Mb(i+1,m) = Mb(i,m) * cosd(FAb);
        %Perform Bloch-McConnell simulation
        [Mb(i+1,m),Mc(i+1,m),Mbtoc(i+1,m)] = bmsim(Mb(i+1,m),Mc(i+1,m),kbc,kcb,TR,T1,np);
    end
    %Calculate pH by summing signal values, correcting for flip angle 
    %differences, and inserting into the Henderson-Hasselbalch equation
    %(iterate for each different value of N)
    %
    %Correct each acquisition for flip angle losses, accounting for exchange
%    Ccorr_factor = (BtoC + 1) * cosd(FAc) - BtoC * pct_exc * (cosd(FAc) - cosd(FAb));
%    Bcorr_factor = (BtoC + 1) * cosd(FAb) - pct_exc * (cosd(FAc) - cosd(FAb));
%    Sc_corr(:,m) = Sc(:,m) ./ (Ccorr_factor .^ [0:size(N,2) - 1]'); 
%    Sb_corr(:,m) = Sb(:,m) ./ (Bcorr_factor .^ [0:size(N,2) - 1]');
    for k=1:size(N,2)
%        Scsum(k,m) = sum(Sc_corr(1:k,m));
%        Sbsum(k,m) = sum(Sb_corr(1:k,m));
        Scsum(k,m) = sum(Sc(1:k,m));
        Sbsum(k,m) = sum(Sb(1:k,m));
        Mbtocsum(k,m) = sum(Mbtoc(1:k,m));
        Sbtocsum(k,m) = sum(Sbtoc(1:k,m));
        %Calculate pH by taking signal ratio and correcting for flip angle
        pH_cal(k,m) = pKa + log10(Sbsum(k,m)/Scsum(k,m)*sind(FAc)/sind(FAb));
    end
    %Calculate SNRs of bicarb and CO2 upon signal averaging
    SNRc(m,:)=Scsum(:,m)'./sqrt(N);
    SNRb(m,:)=Sbsum(:,m)'./sqrt(N);   
    %Determine the value of N that gives max SNR for both CO2 and bicarb
    Nmax(m,1)=find(SNRc(m,:)==max(SNRc(m,:)));
    Nmax(m,2)=find(SNRb(m,:)==max(SNRb(m,:)));
    %Normalize signal to 1st excitation (per pH value), store in separate
    %vector
    Sb_norm(:,m) = Sb(:,m) / max(Sb(:,m));
    Sc_norm(:,m) = Sc(:,m) / max(Sc(:,m));
end

%Generate legend labels
%
for n=1:size(pH,2)
    leg(n,:)=['pH ' num2str(pH(n),'%.1f')];
end

%Calculate signal ratio between first and final acquisitions (averaged over
%all pH values)
%
Sc_ratio_low = Sc(end,1)/Sc(1,1);
Sc_ratio_high = Sc(end,end)/Sc(1,end);
Sb_ratio_low = Sb(end,1)/Sb(1,1);
Sb_ratio_high = Sb(end,end)/Sb(1,end);

%Print results
%
disp(['For BiC ' num2str(FAb,'%2.2f') ' degree flip, CO2 ' num2str(FAc,'%2.2f') ' degree flip, ' num2str(Nexc) ' excitations:'])
%disp(['pH ' num2str(pH(1),2) ': Starting CO2 signal is ' num2str(Sc(1,1),3) ', starting BiC signal is ' num2str(Sb(1,1),3) '\n'])
%disp(['pH ' num2str(pH(end),2) ': Starting CO2 signal is ' num2str(Sc(1,end),3) ', starting BiC signal is ' num2str(Sb(1,end),3) ' \n'])
%disp(['pH ' num2str(pH(1),2) ': Final CO2 signal is ' num2str(Sc(end,1),3) ', final BiC signal is ' num2str(Sb(end,1),3) '\n'])
%disp(['pH ' num2str(pH(end),2) ': Final CO2 signal is ' num2str(Sc(end,end),3) ', final BiC signal is ' num2str(Sb(end,end),3) ' \n'])
disp(['Starting CO2 signal (pH ' num2str(pH(1),2) ') is ' num2str(Sc(1,1),3)])
disp(['Starting CO2 signal (pH ' num2str(pH(end),2) ') is ' num2str(Sc(1,end),3)])
disp(['Starting BiC signal (pH ' num2str(pH(1),2) ') is ' num2str(Sb(1,1),3)])
disp(['Starting BiC signal (pH ' num2str(pH(end),2) ') is ' num2str(Sb(1,end),3)])
disp(['Final CO2 signal (pH ' num2str(pH(1),2) ') is ' num2str(Sc(end,1),3)])
disp(['Final CO2 signal (pH ' num2str(pH(end),2) ') is ' num2str(Sc(end,end),3)])
disp(['Final BiC signal (pH ' num2str(pH(1),2) ') is ' num2str(Sb(end,1),3)])
disp(['Final BiC signal (pH ' num2str(pH(end),2) ') is ' num2str(Sb(end,end),3)])
disp(['Final/starting signal (pH ' num2str(pH(1),2) ') is ' num2str(Sc_ratio_low,3) ' (same for both BiC and CO2)'])
disp(['Final/starting signal (pH ' num2str(pH(end),2) ') is ' num2str(Sc_ratio_high,3) ' (same for both BiC and CO2)'])

%Plot signal per excitation
%
figure
subplot(1,2,1)
plot(N,Sc_norm')
title(['CO_2 Signal per Acquisition, BiC ' num2str(FAb,'%2.2f') '{\circ} flip, CO_2 ' num2str(FAc,'%2.2f') '{\circ} flip'])
xlabel('Acquisition #')
ylabel('Signal')
legend(leg)
subplot(1,2,2)
plot(N,Sb_norm')
title(['BiC Signal per Acquisition, BiC ' num2str(FAb,'%2.2f') '{\circ} flip, CO_2 ' num2str(FAc,'%2.2f') '{\circ} flip'])
xlabel('Acquisition #')
ylabel('Signal')
legend(leg)

%Plot calculated pH vs actual pH (from average BiC and CO2 signals over all
%excitations)
%
figure
plot(pH,pH_cal(end,:))
title(['Calculated pH vs Actual pH, BiC ' num2str(FAb,'%2.2f') '{\circ} flip, CO_2 ' num2str(FAc,'%2.2f') '{\circ} flip'])
xlabel('Actual pH')
ylabel('Calculated pH')
%figure
%subplot(1,2,1)
%plot(N,SNRc)
%title('Signal-Averaged CO_2 SNR vs # of Acquisitions')
%xlabel('Number of Acquisitions')
%ylabel('SNR')
%legend(leg)
%subplot(1,2,2)
%plot(N,SNRb)
%title('Signal-Averaged BiC SNR vs # of Acquisitions')
%xlabel('Number of Acquisitions')
%ylabel('SNR')
%legend(leg)
%subplot(1,3,3)
%plot(pH,Nmax)
%title('Optimal Signal Averaging vs pH')
%xlabel('pH')
%ylabel('Number of Acquisitions for Max SNR')
%legend('CO_2','BiC',2)

% %Plot net %'s of magnetization (bicarb, CO2, total) that transferred from 
% % bicarb over to CO2 over all excitations, as a function of pH
% %
% figure; hold on
% plot(pH,Mbtocsum(end,:) ./ Mb(1,:) * 100);
% plot(pH,Mbtocsum(end,:) ./ Mc(1,:) * 100);
% plot(pH,Mbtocsum(end,:) * 100);
% hold off
% title(['% Bicarb -> CO2 Magnetization Transfer, BiC ' num2str(FAb,'%2.2f') ...
%     '{\circ} flip, CO_2 ' num2str(FAc,'%2.2f') '{\circ} flip'])
% xlabel('pH')
% ylabel('% of Initial Magnetization')
% legend('% of Initial Bicarb','% of Initial CO2','% of Initial Total')

% %Plot net %'s of total CO2 signal arising from bicarb->CO2 exchange over
% %all excitations, as a function of pH
% %
% figure; hold on
% plot(pH,Sbtocsum(end,:) ./ Sbsum(end,:) * 100);
% plot(pH,Sbtocsum(end,:) ./ Scsum(end,:) * 100);
% plot(pH,Sbtocsum(end,:) ./ (Sbsum(end,:) + Scsum(end,:)) * 100);
% hold off
% title(['% CO2 Signal Arising from BiC->CO2 Exchange, BiC ' num2str(FAb,'%2.2f') ...
%     '{\circ} flip, CO_2 ' num2str(FAc,'%2.2f') '{\circ} flip'])
% xlabel('pH')
% ylabel('% of Signal')
% legend('% of Bicarb','% of CO2','% of Total')

%Plot factor increase of CO2 signal over what it would be if there were no
%exchange
%
for m=1:size(pH,2)
    Scsum_noex(m) = sum(cosd(FAc) .^ (N-1) * Mc(1,m) * sind(FAc));
end
ex_enhfac = Scsum(end,:) ./ Scsum_noex;
figure;
plot(pH,ex_enhfac);
title(['% CO2 Signal Enhancement Over No Exchange, BiC ' num2str(FAb,'%2.2f') ...
    '{\circ} flip, CO_2 ' num2str(FAc,'%2.2f') '{\circ} flip'])
xlabel('pH')
ylabel('Enhancement Factor')

%Simulate PSFs due to sampling (Note: if exchange is fast enough, BiC and 
% CO2 are broadened equally for a given pH)
[bpsf,bFWHM] = BiCpsf(Sb,pH,FAb,FAc,1);
[cpsf,cFWHM] = BiCpsf(Sc,pH,FAb,FAc,0);
end

function [Maf,Mbf,Matob] = bmsim(Ma0,Mb0,kab,kba,TR,T1,np)
M = zeros (2 , np);
M(:,1) = [Ma0; Mb0];
R1a = 1 / T1;
R1b = R1a;
dtp = TR / (np-1);
Matob = 0; %magnetization that goes from metabolite A to metabolite B
for i = 1 : np-1
     K = [-1*(R1a+kab)  kba             ;...
          kab           -1*(R1b+kba)    ];
     M(:,i+1) = expm(K * dtp) * M(:,i); %exponential model
     Matob = Matob + M(1,i) * (exp(kab * dtp) - 1) - M(2,i) * (exp(kba * dtp) - 1); 
        %net magnetization that went from A -> B over dtp
end
Maf = M(1,end);
Mbf = M(2,end);
end