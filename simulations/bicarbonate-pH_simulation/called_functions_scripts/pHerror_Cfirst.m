% pHerror_Cfirst.m
% 
% Script that determines the error in calculated pH as a function of flip 
% angle. Sequential excitation is performed on CO2, then bicarbonate (BiC).
% Reconstruction assumes no magnetization exchange between BiC and CO2 
% inbetween excitations. Acquisition is single-shot (1 excitation per 
% resonance).
%
clear all
% Initialize variables
%
FAc = 0 : 5 : 90; % array of flip angles on CO2 (degrees)
FAb = 90; % array of flip angles on BiC (degrees)
pH = 6.5 : 0.1 : 7.4; % pH values used in simulation
ex_rate = 10000000; % exchange rate between BiC and CO2 (1/s)

TR = .154; % repetition time between excitations (s)
T1 = 10; % T1 of BiC and CO2 (s) (assumed equal)
pKa = 6.17; % pKa of BiC-CO2 in vivo at 37 degC
pct_exc = 1 - exp(-1 * TR * ex_rate); % percent exchange completion

Sb = zeros(length(FAc),length(pH));
Sc = Sb;
Mb = Sb;
Mc = Sb;
pH_cal = Sb;
pH_error = Sb;

% Initialize outermost loop, which will iterate each pH value
%
for m = 1 : length(pH)
    % Calculate ratio of bicarb:CO2
    BtoC = 10 ^ (pH(m) - pKa);
    % Initial values of Mb and Mc are calculated using the
    % specified pH and the Henderson-Hasselbalch equation. Total equals 1
    Mb_0 = 1 / (1 / BtoC + 1);
    Mc_0 = 1 / BtoC / (1 / BtoC + 1);
    % Initialize loop, which will calculate the z-magnetization and the signal
    % magnitudes obtained for each flip angle
    for i = 1 : length(FAc)
        % CO2 pulse: obtain signal (xy) magnetization for CO2, adjust 
        % z-magnetization
        Sc(i,m) = Mc_0 * sind(FAc(i));
        Mc(i,m) = Mc_0 * cosd(FAc(i));
        % Magnetization exchange: redistribute total z-magnetization 
        % according to pH
        Mtot = Mc(i,m) + Mb_0;
        Mc_inf = Mtot * 1 / BtoC / (1 / BtoC + 1);
        Mb_inf = Mtot * 1 / (1 / BtoC + 1);
        Mc(i,m) = Mc(i,m) + (Mc_inf - Mc(i,m)) * pct_exc;
        Mb(i,m) = Mb_0 + (Mb_inf - Mb_0) * pct_exc;
        % T1 loss: adjust magnetizations
        Mb(i,m) = Mb(i,m) * exp(-1 * TR / T1);
        Mc(i,m) = Mc(i,m) * exp(-1 * TR / T1);
        % BiC pulse: obtain signal (xy) magnetization for BiC, adjust 
        % z-magnetization
        Sb(i,m) = Mb(i,m) * sind(FAb);
        Mb(i,m) = Mb(i,m) * cosd(FAb);
    end
    % Calculate pH by correcting for flip angle differences and inserting 
    % into the Henderson-Hasselbalch equation (iterate for each different 
    % value of N). Then calculate error from true pH
    pH_cal(:,m) = pKa + log10(squeeze(Sb(:,m)) ./ squeeze(Sc(:,m)) .* sind(FAc') ./ sind(FAb));
    pH_error(:,m) = pH_cal(:,m) - pH(m);
end

% Generate legend labels
%
for n = 1 : size(pH,2)
    leg(n,:) = ['pH ' num2str(pH(n),'%.1f')];
end

% Plot (calculated pH - actual pH) as a function of CO2 flip angle
%
%figure
plot(FAc , pH_error, '--')
%title('pH Error vs CO_2 Flip Angle')
%xlabel('CO_2 Flip Angle (degrees)')
%ylabel('Calculated pH - True pH')
%legend(leg)

% Plot CO2 signal as a function of CO2 flip angle
%
figure
plot(FAc , (Sc ./ Sb)) 
title('CO_2:BiC Signal Ratio vs CO_2 Tip Angle (CO_2 first, BiC tip angle = 90°)')
xlabel('CO_2 Flip Angle (degrees)')
ylabel('CO_2 Signal (arb. units)')
legend(leg)