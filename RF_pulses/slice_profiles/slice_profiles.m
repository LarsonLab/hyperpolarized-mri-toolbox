function [cf , s, gcf, out] = slice_profiles(rf_in, flips, g, verbose)
% This function calculates the correction factors to
% retrospectively apply to HP data or prospectively apply to the slice-select
% gradient to account for slice-profile effects.
% This assumes RF utilization is the main source of HP decay (i.e. TR < T1)
% See Deppe et al 2010 for more details.
%
% INPUTS
% rf_in: This is the RF waveform played out during excitation. The gradient 
% is assumed to be constant during excitation unless explicitly given.
%
% flips: Flip schedule (in radians)
%
% g: gradient waveform (OPTIONAL). If excluded, defaults to a 1 G/cm trapezoid.
%
% verbose: (OPTIONAL). Display diagnostic/final figures, set to 1 to display.
%
% OUTPUTS:
% cf: The correction factor for the HP data. Data must be divided by this
% value to correct for slice-profile effects.
%
% s: Desired signal response based on the given flip schedule.
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


%Initialize variables for simulation
Nrf = length(rf_in);
Nisodelay = Nrf/2;
N = length(flips);
Nsim = 2000;
y = linspace(-1/2, 1/2, Nsim);

if nargin < 3
    g = ones(Nrf,1);
    verbose = 0;
end

if isempty(g)
    g = ones(Nrf,1);
end

if nargin < 4
    verbose = 0;
end

if max(flips) > pi/2
    warning('max(flip) > pi/2, need flips in radians, not degrees')
    beep
    pause(3)
end

%% 1.) Simulate effects for a constant gradient, constant RF
% Estimate signal for all N pulses, accounting for RF decay
s(1) = sin(flips(1));
for ii = 2:N
    s(ii) = sin(flips(ii)).*prod(cos(flips(1:ii-1)));
end
s = s./s(1);

%Scale the RF pulse by the flip angle
for ii = 1:N
    rf_spsp1(:,ii) = rf_in./sum(rf_in).*flips(ii);
end

%First RF pulse
scale = ones(1,N); %same gradient for all
[a b] = abr([rf_spsp1(1:Nrf,1);0], [scale(1).*g;-Nisodelay], y);
Mxyslr_constant(1:Nsim,1) = 2*conj(a).*b;
Mzslr_constant(1:Nsim,1) = 1 - 2*conj(b).*b;

%Desired profile is the first pulse (free of HP effects) that decays by the
%signal profile, calculated above
Mxyslr_desired = Mxyslr_constant;
for ii = 2:N
    Mxyslr_desired(:,ii) = Mxyslr_desired(:,1).*s(ii);
end

%Now the remaining pulses, which will show effects of HP magnetization
for n = 2:N-1
    [a b] = abr([rf_spsp1(1:Nrf,n);0], [scale(n).*g;-Nisodelay], y);
    Mxyslr_constant(1:Nsim,n) = 2*conj(a).*b.* Mzslr_constant(1:Nsim, n-1);
    Mzslr_constant(1:Nsim,n) = ( 1 - 2*conj(b).*b) .* Mzslr_constant(1:Nsim,n-1);
end
for n = N
    [a b] = abr([rf_spsp1(1:Nrf,n);0],[scale(n).*g;-Nisodelay],  y);
    Mxyslr_constant(1:Nsim,n) = 2*conj(a).*b .* Mzslr_constant(1:Nsim, n-1);
    Mzslr_constant(1:Nsim,n) = ( 1 - 2*conj(b).*b ) .* Mzslr_constant(1:Nsim,n-1);
end

%Calculate signal correction factor here
s_spsp = sum(abs(Mxyslr_constant))./max(sum(abs(Mxyslr_constant)));
s_spsp = s_spsp./s_spsp(1); %normalize to 1st value
cf = s_spsp./s;
ratio_uncorr = sum((abs(Mxyslr_constant(1:Nsim,N))))./sum(abs(Mxyslr_constant(1:Nsim,N-1)));

%% 2.) With a constant RF, update the gradients to provide desired response
%%%Calculate gradient correction factor (gcf) here, to yield desired
%%%response
iter = 0;
scale = 1*ones(1,N); %same gradient for all
scale_prior = 1*ones(1,N); %same gradient for all
error = 1;
min_error = Inf;

[a b] = abr([rf_spsp1(1:Nrf,1);0], [scale(1).*g;-Nisodelay], y);
Mxyslr_sim_grad_update(1:Nsim,1) = 2*conj(a).*b;
Mzslr_sim_grad_update(1:Nsim,1) = 1 - 2*conj(b).*b;

%fixed point iteration; break after 20 iterations or less than 1% change
while iter < 20 && error > 1e-2;
    
    for n = 2:N-1
        [a b] = abr([rf_spsp1(1:Nrf,n);0], [scale(n).*g;-Nisodelay], y);
        Mxyslr_sim_grad_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_grad_update(1:Nsim, n-1);
        Mzslr_sim_grad_update(1:Nsim,n) = ( 1 - 2*conj(b).*b) .* Mzslr_sim_grad_update(1:Nsim,n-1);
    end
    for n = N
        [a b] = abr([rf_spsp1(1:Nrf,n);0],[scale(n).*g;-Nisodelay],  y);
        Mxyslr_sim_grad_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_grad_update(1:Nsim, n-1);
        Mzslr_sim_grad_update(1:Nsim,n) = ( 1 - 2*conj(b).*b ) .* Mzslr_sim_grad_update(1:Nsim,n-1);
    end
    
    %Calc AUC ratio for 90:45deg, should be 1 if correctly accounted for
    clc
    f = abs(Mxyslr_sim_grad_update(1:Nsim,N));
    g1 = abs(Mxyslr_sim_grad_update(1:Nsim,N-1));
    dx = 0.001;
    area1 = sum(f.*dx);
    area2 = sum(g1.*dx);
    ratio_corr_grad = abs(area1/area2);
    
    scale_prior = scale;
    area = sum(abs(Mxyslr_sim_grad_update.*dx),1);
    scale = area./area(end);
    scale = (area./max(area(:,1)))./s;
    scale = scale.*scale_prior;
    
    if iter == 0
        error_initial = norm(scale-scale_prior,2);
    end
    
    iter = iter + 1;
    clc
    
    if verbose == 1
        figure(1234)
        drawnow
        subplot(121)
        hold off
        plot(1:n,sum(abs(Mxyslr_sim_grad_update),1)./max(sum(abs(Mxyslr_sim_grad_update(:,1)),1)))
        hold on
        plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired),1)),'--r')
        axis square
        xlim([1 n])
        ylim([0 1.5])
        xlabel('RF #')
        ylabel('Relative Signal')
        legend('Desired','Simulated')
        
        subplot(122)
        plot(scale,'-o')
        axis square
        xlabel('RF #')
        ylabel('Gradient Scaling Factor')
        xlim([1 n])
        
        drawnow
    end

    error_prior = error;
    error = norm(scale-scale_prior,2)./error_initial
    if error < error_prior
        min_error = error;
    end
    if iter > 1 && error > 1.1.*min_error %short circuit
        display('short circuit')
        error = 0;
    end
end

f = abs(Mxyslr_sim_grad_update(1:Nsim,N));
g = abs(Mxyslr_sim_grad_update(1:Nsim,N-1));
dx = 0.001;
area1 = sum(f.*dx);
area2 = sum(g.*dx);
ratio_corr_grad = abs(area1/area2);
gcf = scale;
%% 3.) Keep the gradients constant, update the RF pulse

%Given initial Mz and desired Mxy, get tip profiles!
clear Mxy Mz theta

Mxy = abs(Mxyslr_desired);
Mz(:,1) = ones(length(Mxy),1);
theta(:,1) = asin(Mxy(:,1));

for ii = 2:N
    Mz(:,ii) = Mz(:,ii-1).*cos(theta(:,ii-1));
    theta(:,ii) = asin(Mxy(:,ii)./Mz(:,ii));
end

%If imaginary, we're trying to deliver unrealistic magnetization
theta = real(theta);

%Will crash MATLAB if give too many points for b2a/ab2rf
theta = downsample(theta,round(length(theta)/Nrf));
Nrf = length(theta);

%%%Use tip profile to generate RF pulse
for n = 1:N
    bs(1:Nrf,n) = ifftc(sin(theta(1:Nrf,n)/2));
    as(1:Nrf,n) = b2a(bs(1:Nrf,n));
    rfs(1:Nrf,n) = ab2rf(as(1:Nrf,n),(bs(1:Nrf,n)));
    
    %Window the pulse here to prevent unrealstic excitation
    rfs(1:Nrf,n) = rfs(1:Nrf,n).*hamming(length(rfs));
end

rfs = real(rfs); %Potential issue with complex pulses, need to look into further

%%% Simulate Mxy/Mz with updated RF pulses
z = linspace(0.3,0.7,Nsim);
grad_scale = 8;
[a b] = abr([rfs(1:Nrf,1);0], [grad_scale*ones(Nrf,1);-Nisodelay], y);
Mxyslr_sim_rf_update(1:Nsim,1) = 2*conj(a).*b;
Mzslr_sim_rf_update(1:Nsim,1) = 1 - 2*conj(b).*b;

for n = 2:N-1
    [a b] = abr([rfs(1:Nrf,n);0], [grad_scale*ones(Nrf,1);-Nisodelay], y);
    Mxyslr_sim_rf_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_rf_update(1:Nsim, n-1);
    Mzslr_sim_rf_update(1:Nsim,n) = ( 1 - 2*conj(b).*b) .* Mzslr_sim_rf_update(1:Nsim,n-1);
end
for n = N
    [a b] = abr([rfs(1:Nrf,n);0],[grad_scale*ones(Nrf,1);-Nisodelay],  y);
    Mxyslr_sim_rf_update(1:Nsim,n) = 2*conj(a).*b  .* Mzslr_sim_rf_update(1:Nsim, n-1);
    Mzslr_sim_rf_update(1:Nsim,n) = ( 1 - 2*conj(b).*b ) .* Mzslr_sim_rf_update(1:Nsim,n-1);
end

ratio_corr_rf = sum((abs(Mxyslr_sim_rf_update(1:Nsim,N))))./sum(abs(Mxyslr_sim_rf_update(1:Nsim,N-1)));

%Fill in structure for output
out.mxy_uncorr = Mxyslr_constant;
out.mxy_desired = Mxyslr_desired;
out.mxy_gradcorr = Mxyslr_sim_grad_update;
out.mxy_rfcorr = Mxyslr_sim_rf_update;
out.rf = rfs;

%uncomment below if you wish to save the workspace
% save mag.mat

%% Here is where we'll plot ALL of the results
if verbose == 1
    
    %determine limits for display
    ll = find(abs(Mxyslr_constant(1:Nsim/2,end)) < 0.005.*max(abs(Mxyslr_constant(:,end))),1,'last');
    ul = find(abs(Mxyslr_constant(Nsim/2:Nsim,end)) < 0.005.*max(abs(Mxyslr_constant(:,end))),1,'first');
    ul = Nsim/2 + ul;
    
    %if testing w/ no gradient (g=0), ensure it still works
    if isempty(ll) || isempty(ul)
        ll = 1;
        ul = Nsim;
    else
        ll = ll-100;
        ul = ul+100;
    end
    
    %Stacked slice profiles
    h = figure;
    subplot(131)
    plot(y, abs(Mxyslr_constant)) %look at magnitude, not imaginary component
    title('Slice Profile: Constant Gradient')
    xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
    % xlim([-0.25 0.25])
    xlim([y(ll) y(ul)])
    
    subplot(132)
    plot(y, abs(Mxyslr_sim_grad_update)) %look at magnitude, not imaginary component
    title('Slice Profile: Gradient Update')
    xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
    % xlim([-0.25 0.25])
    xlim([y(ll) y(ul)])
    
    subplot(133)
    plot(y, abs(Mxyslr_sim_rf_update)) %look at magnitude, not imaginary component
    title('Slice Profile: RF Update')
    xlabel('Normalized Frequency'), ylabel('Signal Magnitude')
    % xlim([-0.25 0.25])
    xlim([y(ll) y(ul)])
    set(h,'position',[0 0 1200 800]);
    
    
    %4x2 plot of results
    zmax = max(abs(Mxyslr_constant(:)));
    h = figure;
    subplot(241)
    surf(abs(Mxyslr_desired));
    shading interp
    title('Desired Response')
    view(-72,30)
    axis square
    xlim([1 n])
    ylim([250 750]+500)
    ylim([ll ul])
    zlim([0 zmax])
    
    subplot(242)
    surf(abs(Mxyslr_constant));
    shading interp
    title('Constant External RF')
    view(-72,30)
    axis square
    xlim([1 n])
    ylim([250 750]+500)
    ylim([ll ul])
    zlim([0 zmax])
    
    subplot(243)
    surf(abs(Mxyslr_sim_grad_update));
    shading interp
    title('Constant RF, Update Gradients')
    view(-72,30)
    axis square
    xlim([1 n])
    ylim([250 750]+500)
    ylim([ll ul])
    zlim([0 zmax])
    
    subplot(244)
    surf(abs(Mxyslr_sim_rf_update));
    shading interp
    title('Update RF')
    view(-72,30)
    axis square
    xlim([1 n])
    ylim([250 750]+500)
    ylim([ll ul])
    zlim([0 zmax])
    
    subplot(245)
    %we're just plotting s here (calc above)
    plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)))
    hold on
    plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
    axis square
    xlim([1 n])
    ylim([0 1.5])
    xlabel('RF #')
    title('Normalized Area')
    
    subplot(246)
    % plot(1:n,sum(abs(Mxyslr_constant),1)./max(sum(abs(Mxyslr_constant),1)))
    plot(1:n,sum(abs(Mxyslr_constant),1)./max(sum(abs(Mxyslr_constant(:,1)),1)))
    hold on
    plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
    % plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired),1)),'--r')
    axis square
    xlim([1 n])
    ylim([0 1.5])
    xlabel('RF #')
    title('Normalized Area')
    
    subplot(247)
    plot(1:n,sum(abs(Mxyslr_sim_grad_update),1)./max(sum(abs(Mxyslr_sim_grad_update(:,1)),1)))
    hold on
    plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
    axis square
    xlim([1 n])
    ylim([0 1.5])
    xlabel('RF #')
    title('Normalized Area')
    
    subplot(248)
    plot(1:n,sum(abs(Mxyslr_sim_rf_update),1)./max(sum(abs(Mxyslr_sim_rf_update(:,1)),1)))
    hold on
    plot(1:n,sum(abs(Mxyslr_desired),1)./max(sum(abs(Mxyslr_desired(:,1)),1)),'--r')
    axis square
    xlim([1 n])
    ylim([0 1.5])
    xlabel('RF #')
    title('Normalized Area')
    set(h,'position',[0 0 1600 1067]);

end

clc
sprintf(' Desired ratio for final two pulses: %1.4f\n Uncorrected ratio: %1.4f\n RF corrected ratio: %1.4f\n Gradient corrected ratio: %1.4f',s(N)./s(N-1),ratio_uncorr,ratio_corr_rf,ratio_corr_grad)
end