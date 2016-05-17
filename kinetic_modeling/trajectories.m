function x2 = trajectories( kPL, x1, Mzscale, R1P, R1L , TR )

N = length(x1); 

x2 = zeros(1, N); 
for t=1:N-1
    x2(t+1) = (-(kPL*exp((- R1P - kPL)*TR) - kPL*exp(-R1L*TR))/(R1P - R1L + kPL))*Mzscale(1, t)*x1(t)  +  exp(-R1L*TR)*Mzscale(2, t)*x2(t); 
end

end

