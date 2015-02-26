% Sample variable flip angle generation code
% 
% Author: Peder E. Z. Larson
%
% (c)2014-2015 The Regents of the University of California.
% All Rights Reserved.

%% example for pyruvate (substrate) and lactate (product) 

N = 12; % number of pulses
TR = 3; % repetition time
kPL = 0.05;  % conversion rate from pyruvate to lactate
R1P = 1/35; R1L = 1/30;  % relaxation rates of pyruvate and lactate
Minit = [1; 0.05];  % initial magnetization of Pyr, Lac

% flips_XXX - first row is pyruvate flip angles, second row is lactate flip
% angles

% maximum lactate SNR scheme 
%  (invariant with pyruvate flip angles, so adjust as needed to get 
%   sufficient pyruvate SNR)
flips_maxsnr = [vfa_const_amp(N, pi/2, exp(-TR * (R1P + kPL))); ...
         vfa_opt_signal(N, exp(-TR * R1L ))];

% T1-effective scheme
flips_t1eff = [vfa_const_amp(N, pi/2, exp(-TR * (R1P + kPL))); ...
         vfa_const_amp(N, pi/2, exp(-TR * (R1L - kPL)))];

% const-signal scheme - optimize for constant signal amplitude
%   with the chosen parameters this is not realizable, but increasing N or
%   the initial lactate magnetization, Minit(2), will result in a constant
%   amplitude solution
flips_constsignal = mvfa_const_amp(N, TR, kPL, [R1P R1L], Minit);

[Mxy_maxsnr, Mz_maxsnr] = simulate_2site_model(Minit, [R1P R1L], [kPL 0], flips_maxsnr, TR);
[Mxy_t1eff, Mz_t1eff] = simulate_2site_model(Minit, [R1P R1L], [kPL 0], flips_t1eff, TR);
[Mxy_constsignal, Mz_constsignal] = simulate_2site_model(Minit, [R1P R1L], [kPL 0], flips_constsignal, TR);

t = [0:N-1] * TR;
figure
subplot(221)
plot(t, flips_maxsnr(1,:)*180/pi,t, flips_t1eff(1,:)*180/pi, t, flips_constsignal(1,:)*180/pi)
title('Pyruvate'), legend('Max SNR', 'T1-effective', 'const-signal')
xlabel('time (s)'), ylabel('flip angle (degrees)')
subplot(222)
plot(t, flips_maxsnr(2,:)*180/pi,t, flips_t1eff(2,:)*180/pi, t, flips_constsignal(2,:)*180/pi)
title('Lactate')
xlabel('time (s)'), ylabel('flip angle (degrees)')
subplot(223)
plot(t, Mxy_maxsnr(1,:),t, Mxy_t1eff(1,:), t, Mxy_constsignal(1,:))
title('Pyruvate'), xlabel('time (s)'), ylabel('Simulated Signal')
subplot(224)
plot(t, Mxy_maxsnr(2,:),t, Mxy_t1eff(2,:), t, Mxy_constsignal(2,:))
title('Lactate')
xlabel('time (s)'), ylabel('Simulated Signal')

%%
validate = 0;

if validate
    flips_maxsnrnew = flips_maxsnr; flips_t1effnew = flips_t1eff; flips_constsignalnew = flips_constsignal;
    load test_data/vfa_example_data
    if max(abs(flips_maxsnrnew - flips_maxsnr)) > eps | max(abs(flips_t1effnew - flips_t1eff)) > eps | max(abs(flips_constsignalnew - flips_constsignal )) > eps
        error('Validation failed for variable flip angle example')
    else
        disp('Validation passed for variable flip angle example')
    end
end
