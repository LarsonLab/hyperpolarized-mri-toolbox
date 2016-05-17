function [ l1 ] = negative_log_likelihood_rician(kPL, x1, x2, Mzscale, R1P_fixed, R1L_fixed, TR, noise_level)
    %FUNCTION NEGATIVE_LOG_LIKELIHOOD_RICIAN Computes log likelihood for 
    %    compartmental model with Rician noise 
    
    N = size(x1,2); 
    
    % compute trajectory of the model with parameter values 
    x2fit = trajectories(kPL, x1, Mzscale, R1P_fixed, R1L_fixed, TR);
    
    % compute negative log likelihood 
    l1 = 0;
    for t = 1:N
        for k = 1
            l1 = l1 - (...
                log(x2(k, t)) - log(noise_level(t)) ...
                - (x2(k, t)^2 + x2fit(k, t)^2)/(2*noise_level(t)) ...
                + x2(k, t)*x2fit(k, t)/noise_level(t) ...
                + log(besseli(0, x2(k, t)*x2fit(k, t)/noise_level(t), 1))...
            ); 
        end
    end
end
