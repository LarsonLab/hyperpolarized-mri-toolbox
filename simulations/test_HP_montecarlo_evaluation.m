clear all
NMC = 250;  % 100 is ok

ratio_limits = [.5 1.5];


% default experiment values
exp.R1 = [1/25 1/25]; exp.kPL = 0.02; exp.std_noise = 0.01;

% parameters to plot/fit: kPL, noise, arrival, duration, T1L, T1P
Nplot1 = 3; Nplot2 = 2;

% default fitting parameters
for  est_R1L = 0
    for fit_input = 0
        
        R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;
        params_fixed.R1P = R1P_est;
        params_est.kPL = kPL_est;
        if est_R1L
            params_est.R1L = R1L_est;
        else
            params_fixed.R1L = R1L_est;
        end
        % allowing for fit of R1L increases variability substantially
        % R1P minimal change
        % constraining R1L to narrow range maybe reasonable compromise
        
        if fit_input
            fitting.fit_fcn = @fit_kPL_withinput;
            Tarrival_est = 0;    Tbolus_est = 12;  % ... perfect estimates ... how do they perform with variability?
            Rinj_est = 0.1; % ??
            params_est.Tarrival = Tarrival_est; params_est.Rinj = Rinj_est; params_est.Tbolus = Tbolus_est;
        else
            fitting.fit_fcn = @fit_kPL;
        end
        
        

                    % 2D dynamic 10/20 flips
                    Tacq = 90; acq.TR = 5; acq.N = Tacq/acq.TR;
                    Npe = 8; Nall = acq.N * Npe;
                    acq.flips(1:2,1:acq.N) = repmat(acos(cos([10*pi/180; 20*pi/180]).^Npe), [1 acq.N]);
                    
%                    acq.flips = repmat([10*pi/180; 40*pi/180], [1 N]);

            
            acq.Tbolus = 12;
            Mz0 = [0,0];
            
            fitting.params_est = params_est; fitting.params_fixed = params_fixed;
            
            HP_montecarlo_evaluation( acq, fitting, exp );

    end
end


