function out = nsa(df, ne)
%This function calculates the number of signal averages (NSA) for an IDEAL
%acquisition with n species and ne echoes. Echo spacing (dTE) is assumed to
%be constant.
%
%SYNTAX:
%out = nsa(df,ne)
%
%INPUT: 
%df is a vector containing the relative frequencies for the species of
%interest. 
%
%ne is the number of echoes to acquire     
%
%OUTPUT:
%out.NSA: Number of signal averages (i.e. noise performance) for each species
%out.CN: Condition number
%out.dte: echo spacings used to calculate condition number and NSA
%
%EXAMPLE:
%df = [620 414 288 0]; %lac pa-h2o ala pa @ 4.7T
%out = nsa(df,7);
% 
%JWG 2013-03-07

df = df(:);
n = length(df);
t0 = 1e-3; %s
dte = (0.01:0.01:8)*1e-3; %s

NSA = zeros(n,length(dte));
CN = zeros(1,length(dte));

for ii = 1:length(dte)
    for jj = 1:n
        te = t0:dte(ii):t0+(ne-1)*dte(ii);
        te = te(:);
        A = exp(1i*2*pi*te*df.'); % te*df.' is matrix multiplication
        inv_AHA = inv(A'*A);
        NSA(jj,ii) = 1/inv_AHA(jj,jj);
        CN(1,ii) = cond(A,2);
    end
end

%plot results
figure()
subplot(2,1,1),plot(dte*1e3,real(NSA.'))
hold on
title([int2str(ne) ' echo acquisition for ' int2str(length(df)) ' species']);
xlabel('\DeltaTE (ms)')
ylabel('NSA')
ylim([0 ne+1])

% now plot condition #
subplot(2,1,2)
plot(dte*1e3,CN)
ylim([1 10])
xlabel('\DeltaTE (ms)')
ylabel('Condition #')

% save results in structure
out.nsa = NSA;
out.cond_number = CN;
out.dte = dte;
end