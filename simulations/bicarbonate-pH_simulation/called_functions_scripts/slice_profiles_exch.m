function [cf , s, out] = slice_profiles_exch(rf_in, flips, AtoB, TR, kex, scale)
% This function calculates the correction factors to
% retrospectively apply to HP data or prospectively apply to the slice-select
% gradient to account for slice-profile effects.
% This assumes RF utilization is the main source of HP decay (i.e. TR < T1)
% See Deppe et al 2010 for more details.
%
% INPUTS
% rf_in: This is the RF waveform played out during excitation, for both metabolites. 
% The gradient is assumed to be constant during excitation unless explicitly given.
%
% flips: Flip schedule (in radians). Variable set up as flips(metabolite,flip #)
%
% g: gradient waveform (OPTIONAL). If excluded, defaults to a 1 G/cm trapezoid.
%
% verbose: (OPTIONAL). Display diagnostic/final figures, set to 1 to display.
%
% OUTPUTS:
% cf: The correction factor for the HP data. Data must be divided by this
% value to correct for slice-profile effects. Variable set up as cf.a and
% cf.b
%
% s: Desired signal response based on the given flip schedule. Variable set
%   up as s(metabolite,flip #)
%   NOTE: s.*cf is the expected signal in the presence of these
%   slice-effects.
%
% gcf: Gradient scaling factor to produce the desired signal response.
%
% out: Structure containing magnetization profiles and 
% updated RF pulses to yield the desired signal response.
%
%
% EXAMPLE:
% Determine the expected signal and correction factors for a flip-schedule
% consisting of 15deg pulses repeated 25 times
% load gauss.mat
% rf_in = rf; %Gaussian RF; can replace w/ msinc for sharper pulse
% flips = 15.*ones(1,25).*pi/180; %Constant 15deg flip for 25 excitations 
% [cf, s, gcf, out] = slice_profiles(rf_in,flips);
%
% 9/25/18: Modified by D. Korenchan to incorporate chemical exchange
% between two metabolites A and B. Exchange occurs within each z-position 
% in slice (diffusion assumed to be negligible)
% 10/15/18: Modified by D. Korenchan to include specified gradient
% correction factor for each pulse
% 12/21/18: 


%Initialize variables for simulation
Nrf = length(rf_in);
Nisodelay = Nrf/2;
N = length(flips);
Nsim = 4000;
y = linspace(-1/2, 1/2, Nsim);
g = ones(Nrf,1);
%New variables for exchange
kab = kex * 1 / (AtoB + 1);
kba = kex * AtoB / (AtoB + 1);
Ntp = 100; %# of points to simulate over TR for exchange
Npick = Nsim/5; %profile points to simulate for exchange; makes sim faster!
names = [{'a'},{'b'}];

if max(max(flips)) > pi/2
    warning('max(flip) > pi/2, need flips in radians, not degrees')
    beep
    pause(3)
end

%% 1.) Simulate effects for a constant gradient, constant RF
% Estimate signal for all N pulses, accounting for RF decay and exchange
% for kk = 1:2
%     for ii = 2:N
%         s(kk,ii) = sin(flips(kk,ii)).*prod(cos(flips(kk,1:ii-1)));
%     end
% end
z(:,1) = [AtoB;1];
for ii = 2:N
    z(:,ii) = z(:,ii-1) .* cos(flips(:,ii-1));
    [z(1,ii),z(2,ii)] = bmsimprof(z(1,ii),z(2,ii),1,Ntp);
end
s = z .* sin(flips);
s(1,:) = s(1,:)./s(1,1);
s(2,:) = s(2,:)./s(2,1);

%Scale the RF pulse by the flip angle
for kk = 1:2
    name = names{kk};
    for ii = 1:N
        rf_spsp1.(name)(:,ii) = rf_in./sum(rf_in).*flips(kk,ii);
    end
end

%First RF pulse
if nargin < 6
    scale = ones(1,N); %same gradient for all
end
tic
for kk = 1:2
    name = names{kk};
    [a.(name), b.(name)] = abr([rf_spsp1.(name)(1:Nrf,1);0], [scale(1).*g;-Nisodelay], y);
    Mxyslr_constant.(name)(1:Nsim,1) = 2*conj(a.(name)).*b.(name);
    Mzslr_constant.(name)(1:Nsim,1) = 1 - 2*conj(b.(name)).*b.(name);
    if kk == 1 %scale magnetization values based on AtoB ratio
        Mxyslr_constant.(name)(1:Nsim,1) = Mxyslr_constant.(name)(1:Nsim,1) ...
            * AtoB / (AtoB + 1);
        Mzslr_constant.(name)(1:Nsim,1) = Mzslr_constant.(name)(1:Nsim,1) ...
            * AtoB / (AtoB + 1);
    else
        Mxyslr_constant.(name)(1:Nsim,1) = Mxyslr_constant.(name)(1:Nsim,1) ...
            * 1 / (AtoB + 1);
        Mzslr_constant.(name)(1:Nsim,1) = Mzslr_constant.(name)(1:Nsim,1) ...
            * 1 / (AtoB + 1);
    end

%     %Desired profile is the first pulse (free of HP effects) that decays by the
%     %signal profile, calculated above
%     Mxyslr_desired.(name) = Mxyslr_constant.(name);
%     for ii = 2:N
%         Mxyslr_desired.(name)(:,ii) = Mxyslr_desired.(name)(:,1).*s(kk,ii);
%     end
end

%Before simulating remaining pulses, perform exchange between z-profiles 
%via Bloch-McConnell
Nstart = round(Nsim / 2 - Npick / 2) + 1;
Nend = round(Nsim / 2 + Npick / 2);
[Mzslr_constant.a(Nstart:Nend,1),Mzslr_constant.b(Nstart:Nend,1)] = ...
    bmsimprof(Mzslr_constant.a(Nstart:Nend,1),Mzslr_constant.b(Nstart:Nend,1),Npick,Ntp);

%Now the remaining pulses, which will show effects of HP magnetization
for n = 2:N-1
    %First, the excitations
    for kk = 1:2 
        name = names{kk};
        [a.(name), b.(name)] = abr([rf_spsp1.(name)(1:Nrf,n);0], [scale(n).*g;-Nisodelay], y);
        Mxyslr_constant.(name)(1:Nsim,n) = 2*conj(a.(name)).*b.(name).* Mzslr_constant.(name)(1:Nsim, n-1);
        Mzslr_constant.(name)(1:Nsim,n) = ( 1 - 2*conj(b.(name)).*b.(name)) .* Mzslr_constant.(name)(1:Nsim,n-1);
    end
    %Then, the exchange
    [Mzslr_constant.a(Nstart:Nend,n),Mzslr_constant.b(Nstart:Nend,n)] = ...
        bmsimprof(Mzslr_constant.a(Nstart:Nend,n),Mzslr_constant.b(Nstart:Nend,n),Npick,Ntp);
end
%Final excitations
for kk = 1:2
    name = names{kk};
    for n = N
        [a.(name), b.(name)] = abr([rf_spsp1.(name)(1:Nrf,n);0],[scale(n).*g;-Nisodelay],  y);
        Mxyslr_constant.(name)(1:Nsim,n) = 2*conj(a.(name)).*b.(name) .* Mzslr_constant.(name)(1:Nsim, n-1);
        Mzslr_constant.(name)(1:Nsim,n) = ( 1 - 2*conj(b.(name)).*b.(name) ) .* Mzslr_constant.(name)(1:Nsim,n-1);
    end
    %Calculate signal correction factor here
    s_spsp.(name) = sum(abs(Mxyslr_constant.(name)))./max(sum(abs(Mxyslr_constant.(name))));
    s_spsp.(name) = s_spsp.(name)./s_spsp.(name)(1); %normalize to 1st value
    cf.(name) = s_spsp.(name)./s(kk,:);
%     ratio_uncorr.(name) = sum((abs(Mxyslr_constant.(name)(1:Nsim,N))))./sum(abs(Mxyslr_constant.(name)(1:Nsim,N-1)));
end
toc

%Determine FWHM of composite, summed-together profile for each metabolite
for kk = 1:2
    name = names{kk};
    sumprof.(name).profile = abs(sum(Mxyslr_constant.(name)(1:Nsim,:),2));
    hm = max(sumprof.(name).profile) / 2;
    srch = (sumprof.(name).profile - hm) ./ sumprof.(name).profile;
    sumprof.(name).fwhm = find(abs(srch) < 1e-1,1,'last') - find(abs(srch) < 1e-1,1,'first');
end

% %Determine FWHM for each metabolite, each excitation
% for kk = 1:2
%     name = names{kk};
%     for ii = 1:N
%         hm = max(abs(Mxyslr_constant.(name)(1:Nsim,ii))) / 2;
%         srch = (abs(Mxyslr_constant.(name)(1:Nsim,ii)) - hm) ./ abs(Mxyslr_constant.(name)(1:Nsim,ii));
%         fwhm.(name)(ii) = find(abs(srch) < 1e-1,1,'last') - find(abs(srch) < 1e-1,1,'first');
%     end
%     %Normalize all FWHM's to that of 1st excitation
%     sumprof.(name).fwhm = sumprof.(name).fwhm / fwhm.(name)(1);
%     fwhm.(name) = fwhm.(name) / fwhm.(name)(1);
% end
% 
% disp(['Bicarb (1st metabolite): final composite profile is ' ...
%     num2str(sumprof.a.fwhm,'%2.2f') ' times wider than original pulse'])
% disp(['CO2 (2nd metabolite): final composite profile is ' ...
%     num2str(sumprof.b.fwhm,'%2.2f') ' times wider than original pulse'])

%Stacked slice profiles
h = figure;
subplot(121)
plot(y, abs(Mxyslr_constant.a)/max(max(abs(Mxyslr_constant.a)))) %look at magnitude, not imaginary component
hold on; plot(y, sumprof.a.profile/max(sumprof.a.profile),'-.'); hold off;
title('BiC Slice Profile: Constant Gradient')
xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
% xlim([-0.25 0.25])
xlim([y(Nstart) y(Nend)])

subplot(122)
plot(y, abs(Mxyslr_constant.b)/max(max(abs(Mxyslr_constant.b)))) %look at magnitude, not imaginary component
hold on; plot(y, sumprof.b.profile/max(sumprof.b.profile),'-.'); hold off;
title('CO2 Slice Profile: Constant Gradient')
xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
% xlim([-0.25 0.25])
xlim([y(Nstart) y(Nend)])

%Stacked z-profiles
hh = figure;
subplot(121)
plot(y, abs(Mzslr_constant.a)) %look at magnitude, not imaginary component
title('BiC Z-Profile: Constant Gradient')
xlabel('Normalized Frequency'), ylabel('Mz Magnitude')
% xlim([-0.25 0.25])
xlim([y(Nstart) y(Nend)])

subplot(122)
plot(y, abs(Mzslr_constant.b)) %look at magnitude, not imaginary component
title('CO2 Z-Profile: Constant Gradient')
xlabel('Normalized Frequency'), ylabel('Mz Magnitude')
% xlim([-0.25 0.25])
xlim([y(Nstart) y(Nend)])

%Composite profile over all excitations
hhh = figure;
leg = [{'Composite profile'},{'1st excitation'}];
subplot(121)
plot(y,[sumprof.a.profile / max(sumprof.a.profile) ...
    abs(Mxyslr_constant.a(:,1)) / max(abs(Mxyslr_constant.a(:,1)))])
title('BiC Composite Profile')
xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
xlim([y(Nstart) y(Nend)])
legend(leg)

subplot(122)
plot(y,[sumprof.b.profile / max(sumprof.b.profile) ...
    abs(Mxyslr_constant.b(:,1)) / max(abs(Mxyslr_constant.b(:,1)))])
title('CO2 Composite Profile')
xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
xlim([y(Nstart) y(Nend)])
legend(leg)

% %FWHM over each excitation
% hhhh = figure;
% subplot(121)
% plot(1:N,fwhm.a)
% title('BiC Slice FWHM over Excitations')
% xlabel('Excitation #'), ylabel('FWHM, normalized to 1st excitation')
% 
% subplot(122)
% plot(1:N,fwhm.b)
% title('CO2 Slice FWHM over Excitations')
% xlabel('Excitation #'), ylabel('FWHM, normalized to 1st excitation')

%Fill in structure for output
out.mxy_uncorr = Mxyslr_constant;
out.mz_uncorr = Mzslr_constant;
out.signal_total = s_spsp;
out.sumprof = sumprof;
% out.fwhm = fwhm;

% %% 2.) With a constant RF, update the gradients to provide desired response
% %%%Calculate gradient correction factor (gcf) here, to yield desired
% %%%response
% iter = 0;
% scale = 1*ones(1,N); %same gradient for all
% scale_prior = 1*ones(1,N); %same gradient for all
% error = 1;
% min_error = Inf;
% 
% [a b] = abr([rf_spsp1(1:Nrf,1);0], [scale(1).*g;-Nisodelay], y);
% Mxyslr_sim_grad_update(1:Nsim,1) = 2*conj(a).*b;
% Mzslr_sim_grad_update(1:Nsim,1) = 1 - 2*conj(b).*b;
% 
% %fixed point iteration; break after 20 iterations or less than 1% change
% while iter < 20 && error > 1e-2;
%     
%     for n = 2:N-1
%         [a b] = abr([rf_spsp1(1:Nrf,n);0], [scale(n).*g;-Nisodelay], y);
%         Mxyslr_sim_grad_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_grad_update(1:Nsim, n-1);
%         Mzslr_sim_grad_update(1:Nsim,n) = ( 1 - 2*conj(b).*b) .* Mzslr_sim_grad_update(1:Nsim,n-1);
%     end
%     for n = N
%         [a b] = abr([rf_spsp1(1:Nrf,n);0],[scale(n).*g;-Nisodelay],  y);
%         Mxyslr_sim_grad_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_grad_update(1:Nsim, n-1);
%         Mzslr_sim_grad_update(1:Nsim,n) = ( 1 - 2*conj(b).*b ) .* Mzslr_sim_grad_update(1:Nsim,n-1);
%     end
%     
%     %Calc AUC ratio for 90:45deg, should be 1 if correctly accounted for
%     clc
%     f = abs(Mxyslr_sim_grad_update(1:Nsim,N));
%     g1 = abs(Mxyslr_sim_grad_update(1:Nsim,N-1));
%     dx = 0.001;
%     area1 = sum(f.*dx);
%     area2 = sum(g1.*dx);
%     ratio_corr_grad = abs(area1/area2);
%     
%     scale_prior = scale;
%     area = sum(abs(Mxyslr_sim_grad_update.*dx),1);
%     scale = area./area(end);
%     scale = (area./max(area(:,1)))./s;
%     scale = scale.*scale_prior;
%     
%     if iter == 0
%         error_initial = norm(scale-scale_prior,2);
%     end
%     
%     iter = iter + 1;
%     clc
%     
%     if verbose == 1
%         figure(1234)
%         drawnow
%         subplot(121)
%         hold off
%         plot(1:n,sum(abs(Mxyslr_sim_grad_update),1)./max(sum(abs(Mxyslr_sim_grad_update(:,1)),1)))
%         hold on
%         plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired),1)),'--r')
%         axis square
%         xlim([1 n])
%         ylim([0 1.5])
%         xlabel('RF #')
%         ylabel('Relative Signal')
%         legend('Desired','Simulated')
%         
%         subplot(122)
%         plot(scale,'-o')
%         axis square
%         xlabel('RF #')
%         ylabel('Gradient Scaling Factor')
%         xlim([1 n])
%         
%         drawnow
%     end
% 
%     error_prior = error;
%     error = norm(scale-scale_prior,2)./error_initial
%     if error < error_prior
%         min_error = error;
%     end
%     if iter > 1 && error > 1.1.*min_error %short circuit
%         display('short circuit')
%         error = 0;
%     end
% end
% 
% f = abs(Mxyslr_sim_grad_update(1:Nsim,N));
% g = abs(Mxyslr_sim_grad_update(1:Nsim,N-1));
% dx = 0.001;
% area1 = sum(f.*dx);
% area2 = sum(g.*dx);
% ratio_corr_grad = abs(area1/area2);
% gcf = scale;
% %% 3.) Keep the gradients constant, update the RF pulse
% 
% %Given initial Mz and desired Mxy, get tip profiles!
% clear Mxy Mz theta
% 
% Mxy = abs(Mxyslr_desired);
% Mz(:,1) = ones(length(Mxy),1);
% theta(:,1) = asin(Mxy(:,1));
% 
% for ii = 2:N
%     Mz(:,ii) = Mz(:,ii-1).*cos(theta(:,ii-1));
%     theta(:,ii) = asin(Mxy(:,ii)./Mz(:,ii));
% end
% 
% %If imaginary, we're trying to deliver unrealistic magnetization
% theta = real(theta);
% 
% %Will crash MATLAB if give too many points for b2a/ab2rf
% theta = downsample(theta,round(length(theta)/Nrf));
% Nrf = length(theta);
% 
% %%%Use tip profile to generate RF pulse
% for n = 1:N
%     bs(1:Nrf,n) = ifftc(sin(theta(1:Nrf,n)/2));
%     as(1:Nrf,n) = b2a(bs(1:Nrf,n));
%     rfs(1:Nrf,n) = ab2rf(as(1:Nrf,n),(bs(1:Nrf,n)));
%     
%     %Window the pulse here to prevent unrealstic excitation
%     rfs(1:Nrf,n) = rfs(1:Nrf,n).*hamming(length(rfs));
% end
% 
% rfs = real(rfs); %Potential issue with complex pulses, need to look into further
% 
% %%% Simulate Mxy/Mz with updated RF pulses
% z = linspace(0.3,0.7,Nsim);
% grad_scale = 8;
% [a b] = abr([rfs(1:Nrf,1);0], [grad_scale*ones(Nrf,1);-Nisodelay], y);
% Mxyslr_sim_rf_update(1:Nsim,1) = 2*conj(a).*b;
% Mzslr_sim_rf_update(1:Nsim,1) = 1 - 2*conj(b).*b;
% 
% for n = 2:N-1
%     [a b] = abr([rfs(1:Nrf,n);0], [grad_scale*ones(Nrf,1);-Nisodelay], y);
%     Mxyslr_sim_rf_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_rf_update(1:Nsim, n-1);
%     Mzslr_sim_rf_update(1:Nsim,n) = ( 1 - 2*conj(b).*b) .* Mzslr_sim_rf_update(1:Nsim,n-1);
% end
% for n = N
%     [a b] = abr([rfs(1:Nrf,n);0],[grad_scale*ones(Nrf,1);-Nisodelay],  y);
%     Mxyslr_sim_rf_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_rf_update(1:Nsim, n-1);
%     Mzslr_sim_rf_update(1:Nsim,n) = ( 1 - 2*conj(b).*b ) .* Mzslr_sim_rf_update(1:Nsim,n-1);
% end
% 
% ratio_corr_rf = sum((abs(Mxyslr_sim_rf_update(1:Nsim,N))))./sum(abs(Mxyslr_sim_rf_update(1:Nsim,N-1)));
% 
% 
%uncomment below if you wish to save the workspace
% save mag.mat
% 
% %% Here is where we'll plot ALL of the results
% if verbose == 1
%     
%     %determine limits for display
%     ll = find(abs(Mxyslr_constant(1:Nsim/2,end)) < 0.005.*max(abs(Mxyslr_constant(:,end))),1,'last');
%     ul = find(abs(Mxyslr_constant(Nsim/2:Nsim,end)) < 0.005.*max(abs(Mxyslr_constant(:,end))),1,'first');
%     ul = Nsim/2 + ul;
%     
%     %if testing w/ no gradient (g=0), ensure it still works
%     if isempty(ll) || isempty(ul)
%         ll = 1;
%         ul = Nsim;
%     else
%         ll = ll-100;
%         ul = ul+100;
%     end
%     
%     %Stacked slice profiles
%     h = figure;
%     subplot(131)
%     plot(y, abs(Mxyslr_constant)) %look at magnitude, not imaginary component
%     title('Slice Profile: Constant Gradient')
%     xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
%     % xlim([-0.25 0.25])
%     xlim([y(ll) y(ul)])
%     
%     subplot(132)
%     plot(y, abs(Mxyslr_sim_grad_update)) %look at magnitude, not imaginary component
%     title('Slice Profile: Gradient Update')
%     xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
%     % xlim([-0.25 0.25])
%     xlim([y(ll) y(ul)])
%     
%     subplot(133)
%     plot(y, abs(Mxyslr_sim_rf_update)) %look at magnitude, not imaginary component
%     title('Slice Profile: RF Update')
%     xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
%     % xlim([-0.25 0.25])
%     xlim([y(ll) y(ul)])
%     set(h,'position',[0 0 1200 800]);
%     
%     
%     %4x2 plot of results
%     zmax = max(abs(Mxyslr_constant(:)));
%     h = figure;
%     subplot(241)
%     surf(abs(Mxyslr_desired));
%     shading interp
%     title('Desired Response')
%     view(-72,30)
%     axis square
%     xlim([1 n])
%     ylim([250 750]+500)
%     ylim([ll ul])
%     zlim([0 zmax])
%     
%     subplot(242)
%     surf(abs(Mxyslr_constant));
%     shading interp
%     title('Constant External RF')
%     view(-72,30)
%     axis square
%     xlim([1 n])
%     ylim([250 750]+500)
%     ylim([ll ul])
%     zlim([0 zmax])
%     
%     subplot(243)
%     surf(abs(Mxyslr_sim_grad_update));
%     shading interp
%     title('Constant RF, Update Gradients')
%     view(-72,30)
%     axis square
%     xlim([1 n])
%     ylim([250 750]+500)
%     ylim([ll ul])
%     zlim([0 zmax])
%     
%     subplot(244)
%     surf(abs(Mxyslr_sim_rf_update));
%     shading interp
%     title('Update RF')
%     view(-72,30)
%     axis square
%     xlim([1 n])
%     ylim([250 750]+500)
%     ylim([ll ul])
%     zlim([0 zmax])
%     
%     subplot(245)
%     %we're just plotting s here (calc above)
%     plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)))
%     hold on
%     plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
%     axis square
%     xlim([1 n])
%     ylim([0 1.5])
%     xlabel('RF #')
%     title('Normalized Area')
%     
%     subplot(246)
%     % plot(1:n,sum(abs(Mxyslr_constant),1)./max(sum(abs(Mxyslr_constant),1)))
%     plot(1:n,sum(abs(Mxyslr_constant),1)./max(sum(abs(Mxyslr_constant(:,1)),1)))
%     hold on
%     plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
%     % plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired),1)),'--r')
%     axis square
%     xlim([1 n])
%     ylim([0 1.5])
%     xlabel('RF #')
%     title('Normalized Area')
%     
%     subplot(247)
%     plot(1:n,sum(abs(Mxyslr_sim_grad_update),1)./max(sum(abs(Mxyslr_sim_grad_update(:,1)),1)))
%     hold on
%     plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
%     axis square
%     xlim([1 n])
%     ylim([0 1.5])
%     xlabel('RF #')
%     title('Normalized Area')
%     
%     subplot(248)
%     plot(1:n,sum(abs(Mxyslr_sim_rf_update),1)./max(sum(abs(Mxyslr_sim_rf_update(:,1)),1)))
%     hold on
%     plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
%     axis square
%     xlim([1 n])
%     ylim([0 1.5])
%     xlabel('RF #')
%     title('Normalized Area')
%     set(h,'position',[0 0 1600 1067]);
% 
% end
% 
% clc
% sprintf(' Desired ratio for final two pulses: %1.4f\n Uncorrected ratio: %1.4f\n RF corrected ratio: %1.4f\n Gradient corrected ratio: %1.4f',s(N)./s(N-1),ratio_uncorr,ratio_corr_rf,ratio_corr_grad)

%% INTERNAL FUNCTIONS
%
%bmsimprof: Takes two z-magnetization profiles in z and performs
%Bloch-McConnell simulation between them, using the specified TR and first-
%order forward/reverse exchange rates. T1 decay and lateral diffusion along
%the z-axis are neglected.
%
function [Maf,Mbf] = bmsimprof(Ma0,Mb0,Nprofpts,np)
M = zeros (Nprofpts*2 , np);
M(:,1) = [Ma0; Mb0];
dtp = TR / (np-1);
K = [diag(ones(1,Nprofpts)*-kab)    diag(ones(1,Nprofpts)*kba)  ;
    diag(ones(1,Nprofpts)*kab)      diag(ones(1,Nprofpts)*-kba) ];
Ke = expm(K * dtp);
for i = 1 : np-1
     M(:,i+1) = Ke * M(:,i); %exponential model
end
Maf = M(1:Nprofpts,end);
Mbf = M(Nprofpts+1:end,end);
end
end