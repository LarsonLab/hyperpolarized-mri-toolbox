function [Sscale, Mzscale,f] = extract_flip_compensation_profile(root_dir, root_fname, Npe, Nt, bw, Nf);
% [Sscale, Mzscale,f] = extract_flip_compensation_profile(root_dir, root_fname, Npe, Nt, bw, Nf);
%
% Determines the effect on the hyperpolarized magnetization of RF pulses
% across a bandwidth.   Computes the loss of MZ (Mzscale) and the effective
% signal scaling (Sscale) that can be included in kinetic modeling as
% described in:
%   Bahrami N, Swisher CL, von Morze C, Vigneron DB, Larson PE. Kinetic and
%   perfusion modeling of hyperpolarized 13C pyruvate and urea in cancer
%   with arbitrary RF flip angles. Quant Imaging Med Surg 2014;4(1):24-32.
%   Doi: 10.3978/j.issn.2223-4292.2014.02.02
%
% INPUTS:
%   root_dir - root directory where pulses are stored
%   root_fname - root file name where pulses are stored
%   Npe - number of phase encodes per image
%   Npe - number of phase encodes per image
%   Nt - number of time points
%   bw, Nf - bandwidth and number of frequency points to compute response
%       over
%
% OUTPUTS:
%   Sscale - effective signal scaling factors
%   Mzscale - effective Mz scaling factors
%   f - frequency vector of computation
%
% (c) 2016-2017 The Regents of the University of California
% All Rights Reserved.
%
% Originally written by: Peder E. Z. Larson


%% get pulse spectral profile at z=0

f = ([0:Nf-1]/Nf - 0.5) * bw;


for n = 1:Nt*Npe
    
    % create ss_read() ?
    rfrho = signaread(sprintf('%s/%s%03d.rho', root_dir, root_fname, n));
    rfpha = signaread(sprintf('%s/%s%03d.pha', root_dir, root_fname, n));
    
    fid = fopen(sprintf('%s/%s%03d.dat', root_dir, root_fname, n), 'r');
    rfdat = fscanf(fid, '%f \t\t %*s\n', Inf);  
    rfpeakB1= rfdat(10) *4257/1071;
%     fprintf(fid,'%10d \t\t #Spectral-Spatial\n', 1);
%             fprintf(fid,'%10d \t\t #res\n', length(rf));
%             fprintf(fid,'%10d \t\t #pw\n',round(length(rf)*SS_TS*1e6));
%             fprintf(fid,'%10.7f \t\t #nom_flip \n',ang*180/pi);
%             fprintf(fid,'%10.7f \t\t #abswidth \n',abswidth);
%             fprintf(fid,'%10.7f \t\t #effwidth \n',effwidth);
%             fprintf(fid,'%10.7f \t\t #area \n',area);
%             fprintf(fid,'%10.7f \t\t #dtycyc \n',dty_cyc);
%             fprintf(fid,'%10.7f \t\t #maxpw \n',max_pw);
%             gamscale = SS_GAMMA/SS_GAMMA_HYDROGEN; % GE assumes max B1 is for application at 1H
%             fprintf(fid,'%10.7f \t\t #max_b1 \n',max_b1 * gamscale); 
%             fprintf(fid,'%10.7f \t\t #max_int_b1_sqr \n',int_b1_sqr* gamscale^2);
%             fprintf(fid,'%10.7f \t\t #max_rms_b1 \n',rms_b1* gamscale^2);
%             fprintf(fid,'%10.3f \t\t #a_gzs \n',maxg);
%             fprintf(fid,'%10.3f \t\t #nom_thk(mm) \n',thk * gamscale * 10);
%             fprintf(fid,'%10d \t\t #isodelay\n',round(isodelay*1e6));

    fclose(fid);

    Nrf = rfdat(2);
    dt = rfdat(3)/rfdat(2);
    
    % remove header in GE formatted RF pulse file
    rf = rfpeakB1 * rfrho([1:Nrf]+32)/32767 .* exp(i*rfpha([1:Nrf]+32) *pi/32767);

    [mxy(1:Nf,n) mz(1:Nf,n)] = ss_response_mxy(zeros(size(rf)),rf,0,f,dt *1e-6,1071);
 
    
end


%%
Mzscale = zeros(Nf, Nt); Sscale = zeros(Nf, Nt);
for m= 1:Nt
    Iflips = [1:Npe] + (m-1)*Npe;
    Mzscale(:,m) = prod(mz(:,Iflips),2);
    for n = 1:Npe
        Sscale(:,m) = Sscale(:,m) + abs(mxy(:,Iflips(n))) .* prod(mz(:,Iflips(1:n-1)),2);
    end
end
Sscale = Sscale/Npe;


%%

fid = fopen([root_dir '/' root_fname '_flip_profile.dat'],'w');

fprintf(fid, 'dat_type: variable_flip_angle\n');
fprintf(fid, 'dat_version: 1\n');

fprintf(fid,'dat_content_name: %s%s\n', root_dir, root_fname);
fprintf(fid,'num_time_pts: %d\n', Nt);
fprintf(fid,'num_phase_encodes: %d\n', Npe);
fprintf(fid,'profile_bandwidth: %f\n', bw);
fprintf(fid,'profile_num_pts: %d\n', Nf);

for t=1:Nt
    fprintf(fid,'time_point: %d\n',t);
    fprintf(fid,'signal_scaling: ',t);
    fprintf(fid,'%f ',Sscale(:,t));
    fprintf(fid, '\n');
    fprintf(fid,'mz_scaling: ',t);
    fprintf(fid,'%f ',Mzscale(:,t));
    fprintf(fid, '\n');
end

fclose(fid);
    