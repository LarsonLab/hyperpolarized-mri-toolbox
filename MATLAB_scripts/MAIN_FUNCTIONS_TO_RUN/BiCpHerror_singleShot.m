% BiCpHerror_singleShot.m
% 
% Script that determines the error in calculated pH as a function of flip 
% angle for sequential, single-shot imaging on BiC and CO2 (ie. 1 
% excitation per resonance). BiC is excited, and exchange between
% z-magnetization pools is allowed to occur before CO2 is excited. The pH
% is then calculated from the obtained signals, correcting for tip angles
% used. This simulation is then repeated with CO2 excited first (under
% script pHerror_Cfirst.m). 
%
clear all
% Initialize variables
%
FAb = 0 : 5 : 90; % array of flip angles on BiC (degrees)
FAc = 90; % array of flip angles on CO2 (degrees)
pH = 6.4 : 0.1 : 7.6; % pH values used in simulation
ex_rate = 10000000; % exchange rate between BiC and CO2 (1/s)

TR = .154; % repetition time between excitations (s)
T1 = 10; % T1 of BiC and CO2 (s) (assumed equal)
pKa = 6.17; % pKa of BiC-CO2 in vivo at 37 degC
pct_exc = 1 - exp(-1 * TR * ex_rate); % percent exchange completion

Sb = zeros(length(FAb),length(pH));
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
    for i = 1 : length(FAb)
        % BiC pulse: obtain signal (xy) magnetization for BiC, adjust 
        % z-magnetization
        Sb(i,m) = Mb_0 * sind(FAb(i));
        Mb(i,m) = Mb_0 * cosd(FAb(i));
        % Magnetization exchange: redistribute total z-magnetization 
        % according to pH
        Mtot = Mb(i,m) + Mc_0;
        Mc_inf = Mtot * 1 / BtoC / (1 / BtoC + 1);
        Mb_inf = Mtot * 1 / (1 / BtoC + 1);
        Mc(i,m) = Mc_0 + (Mc_inf - Mc_0) * pct_exc;
        Mb(i,m) = Mb(i,m) + (Mb_inf - Mb(i,m)) * pct_exc;
        % T1 loss: adjust magnetizations
        Mb(i,m) = Mb(i,m) * exp(-1 * TR / T1);
        Mc(i,m) = Mc(i,m) * exp(-1 * TR / T1);
        % CO2 pulse: obtain signal (xy) magnetization for CO2, adjust 
        % z-magnetization
        Sc(i,m) = Mc(i,m) * sind(FAc);
        Mc(i,m) = Mc(i,m) * cosd(FAc);
    end
    % Calculate pH by correcting for flip angle differences and inserting 
    % into the Henderson-Hasselbalch equation (iterate for each different 
    % value of N). Then calculate error from true pH
    pH_cal(:,m) = pKa + log10(squeeze(Sb(:,m)) ./ squeeze(Sc(:,m)) .* sind(FAc) ./ sind(FAb'));
    pH_error(:,m) = pH_cal(:,m) - pH(m);
end

% Generate legend labels
%
for n = 1 : size(pH,2)
    leg(n,:) = ['pH ' num2str(pH(n),'%.1f')];
end

% Plot CO2 signal as a function of BiC flip angle
%
figure
plot(FAb , (Sc ./ Sb)) 
title('CO_2:BiC Signal Ratio vs BiC Tip Angle (BiC first, CO_2 tip angle = 90°)')
xlabel('BiC Flip Angle (degrees)')
ylabel('CO_2 Signal (arb. units)')
legend(leg)

% Plot (calculated pH - actual pH) as a function of CO2 flip angle
%
figure
plot(FAb , pH_error)
title('pH Error vs BiC Tip Angle')
xlabel('BiC Flip Angle (degrees)')
ylabel('Calculated pH - True pH')
legend(leg)

hold on
pHerror_Cfirst