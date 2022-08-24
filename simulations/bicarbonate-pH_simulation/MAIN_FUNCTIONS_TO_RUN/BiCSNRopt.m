function BiCSNRopt(infilename)
%BiCSNRopt.m
%
%Used to identify imaging parameters for multi-excitation pH imaging that 
%will maximize imaging SNR while constraining pH error within a specified
%tolerance. The script first performs Bloch-McConnell simulation over
%specified true pH values and parameter combinations, including bicarb/CO2
%tip angles, TR between excitations, and total number of excitations 
%(simulation is skipped if input file name is specified as input). The 
%function then calculates both pH error and SNR values as signal is
%averaged over excitations (assuming Gaussian-distributed, random-phase 
%noise). The function then plots an interactive figure in which CO2 SNR as
%a function of 2 selected parameter values is displayed as a 3D plot, and
%parameter combinations maximizing the minimum BiC/CO2 SNR for each pH
%value while enforcing the pH error cutoff.
%
%INPUTS:
%
%infilename:    String matching the name of the .mat file containing
%               previously-saved simulation data in directory specified by 
%               save_dir variable; if specified, script will load variables 
%               and display interactive figure only
%
%MODIFICATIONS:
%11/28/18:  Adjusted script so that SNR masking would only look at min BiC
%           and min CO2 SNRs over pH range, once both are maximized for Nexc
%11/29/18:  Added display of pH error for each simulated max
%3/25/19:   Changed to specify timestep between points, rather than # of 
%           points, for Bloch-McConnell simulation


%% VALUES TO ADJUST PRIOR TO RUNNING
%
% save_dir = '/home/dkorenchan/Dave_K/MATLAB/Simulation/Saved_data';
save_dir = '/Users/sf865719/Lab/MATLAB/Simulation/HP_BiC/Saved_data';
    %name of directory for loading previously-saved simulation data
npBM = 0;
dtBM = .01; %time step between points for Bloch-McConnell (s)
        
% Variables used in simulation: "givens"
%
pH = [6.4 7.6]; %pH range (min, max) for simulation - phantom
% pH = [6.8 7.6]; %pH range (min, max) for simulation - in vivo
% pH = [7.2 7.0];
% T1 = 10; %T1 of BiC and CO2 (s, assumed equal) - in vivo
T1 = 25; %T1 of BiC and CO2 (s, assumed equal) - in solution, 14 T
% kex = 0;
% kex = 1.56; %exchange rate between BiC and CO2 (1/s) - measured in TRAMP tumors
% kex = .168; %exchange rate between BiC and CO2 (1/s) - in free solution, no CA
% kex = 5.51; %exchange rate between BiC and CO2 (1/s) - with CA at 7.55 ug/mL, pH 6.7
% kex = 4.26; %exchange rate between BiC and CO2 (1/s) - with CA at 7.55 ug/mL, pH 7.8 (calculated)
kex = 4.32; %exchange rate between BiC and CO2 (1/s) - with CA at 7.55 ug/mL, pH 7.6 (calculated)

% Variables used in simulation: sequence parameters to optimize
%
pars.FAb = 1:10; %tip angle, BiC (°)
pars.FAc = 46:2:82; %tip angle, CO2 (°)
pars.TR = 0.1:0.1:2; %repetition time between excitations (s)
pars.Nexc = 1:64; %# of excitations averaged together
    %NOTE: Nexc MUST start @ 1 and go up by 1's! Otherwise won't work


%% CHECK FOR PREVIOUS FILE
%
% If filename specified as input, load and skip simulation
%
clc
switch(nargin)
    case 1
        disp('Input filename detected. Loading .mat file...')
        load([infilename '.mat']);
    otherwise %perform simulation
        disp('No input specified. Performing simulation...')
        %% DEFINE VARIABLES
        %
        home = pwd;
        tab = sprintf('\t');

        disp('Simulation parameters:')
        disp(['minpH' tab 'maxpH' tab 'T1' tab 'k_ex'])
        disp([num2str(pH(1),'%1.1f') tab num2str(pH(end),'%1.1f') tab ...
        num2str(T1,'%2.0f') tab num2str(kex,'%2.2f')])

        if pars.Nexc(1) ~= 1 || pars.Nexc(end) ~= length(pars.Nexc) %catch improper Nexc vector
            error('Nexc not properly set! Must start at 1, increment by 1 each time')
        end

        % Initialize variables (necessary for simulation)
        %
        pKa = 6.17; %pKa of BiC-CO2 in vivo at 37 degC
        BtoC = 10 .^ (pH - pKa); %ratio of bicarb to CO2, equilibrium
        n.FAb = length(pars.FAb);
        n.FAc = length(pars.FAc);
        n.TR = length(pars.TR);
        n.Nexc = length(pars.Nexc);
        n.pH = length(pH);

        %Matrices for simulation
        Mb = zeros(n.FAb,n.FAc,n.TR,n.Nexc,n.pH); %bicarb z-magnetization just prior to sampling
        Mc = Mb; %CO2 z-magnetization just prior to sampling
        Sb = Mb; %summed bicarb signal over all previous excitations
        Sc = Mb; %summed CO2 signal over all previous excitations
        noise = Mb; %noise stdev
        noise0 = 0.0072; %initial noise
        pHcalc = Mb; %calculated pH from signal
        pHerr = Mb; %calculated pH error from summed signals

        % Compute noise matrix
        for ll = 1:n.Nexc
            noise(:,:,:,ll,:) = noise0 * sqrt(ll);
        end


        %% SIMULATION
        %
        f = waitbar(0,'Simulating...','WindowStyle','docked');
        tic
        for mm = 1:n.pH %loop for min/max pH's
            kbc = kex * 1/(BtoC(mm)+1); %bicarb-to-CO2 first-order rate, 1/s 
            kcb = kex * BtoC(mm)/(BtoC(mm)+1); %CO2-to-bicarb first-order rate, 1/s
            Mb(:,:,:,1,mm) = BtoC(mm)/(BtoC(mm)+1);
            Mc(:,:,:,1,mm) = 1/(BtoC(mm)+1);
            for ii = 1:n.FAb %bicarb tip angle loop
                for jj = 1:n.FAc %CO2 tip angle loop
                    for kk = 1:n.TR %TR loop
                        %1st excitation
                        Sb(ii,jj,kk,1,mm) = Mb(ii,jj,kk,1,mm) * sind(pars.FAb(ii));
                        Sc(ii,jj,kk,1,mm) = Mc(ii,jj,kk,1,mm) * sind(pars.FAc(jj));
                        Mb(ii,jj,kk,2,mm) = Mb(ii,jj,kk,1,mm) * cosd(pars.FAb(ii));          
                        Mc(ii,jj,kk,2,mm) = Mc(ii,jj,kk,1,mm) * cosd(pars.FAc(jj));
                        for ll = 2:n.Nexc %Nexc loop
                            %Bloch-McConnell on z-magnetization
                            [Mb(ii,jj,kk,ll,mm),Mc(ii,jj,kk,ll,mm)] = ...
                                bmsim(Mb(ii,jj,kk,ll,mm),Mc(ii,jj,kk,ll,mm),kbc,kcb,pars.TR(kk),T1,dtBM);
                            %'ll'th excitation
                            Sb(ii,jj,kk,ll,mm) = Sb(ii,jj,kk,ll-1,mm) + Mb(ii,jj,kk,ll,mm) * sind(pars.FAb(ii));
                            Sc(ii,jj,kk,ll,mm) = Sc(ii,jj,kk,ll-1,mm) + Mc(ii,jj,kk,ll,mm) * sind(pars.FAc(jj));
                            if ll ~= n.Nexc %don't do on final excitation
                                Mb(ii,jj,kk,ll+1,mm) = Mb(ii,jj,kk,ll,mm) * cosd(pars.FAb(ii));
                                Mc(ii,jj,kk,ll+1,mm) = Mc(ii,jj,kk,ll,mm) * cosd(pars.FAc(jj));
                            end                 
                        end
                    end
                    %Calculate pH, correcting for tip angles
                    pHcalc(ii,jj,:,:,mm) = pKa + log10(Sb(ii,jj,:,:,mm) ./ Sc(ii,jj,:,:,mm) ...
                        * sind(pars.FAc(jj)) / sind(pars.FAb(ii)));
                    waitbar(1/(n.FAb*n.FAc*n.pH)*(n.FAc*n.FAb*(mm-1)+n.FAc*(ii-1)+jj),f,'Simulating...')
                end
            end
        end
        toc
        waitbar(1,f,'Complete!');
        close(f)

        % Save matrices
        %
        clear f %to reduce file size
        basefilename = ['SIMpH' num2str(pH(1),'%1.1f') 'to' num2str(pH(end),'%1.1f') ...
            'Tone' num2str(T1,'%2.0f') 'kex' num2str(kex,'%2.2f')];
        cd(save_dir)
        ctr = 1;
        filename = basefilename;
        while 1 %this loop prevents saving over previous data
            if exist([filename '.mat'],'file')
                ctr = ctr + 1;
                filename = [basefilename '_' num2str(ctr,'%i')];
            else
                break;
            end
        end
        save([filename '.mat'])
        disp(['Workspace saved as ' filename '.mat'])
        cd(home)
end


%% SNR CALCULATION, pH ACCURACY CHECKING
%
% Calculate SNR, pH error matrices
%
SNRb = Sb ./ noise;
SNRc = Sc ./ noise;
% SNRbtoc = SNRb ./ SNRc; %BiC/CO2 SNR ratio for each set of parameters, given pH
SNRbtocchk = repmat(min(max(SNRb,[],4),[],5),1,1,1,n.Nexc,n.pH) ./ ...
    repmat(min(max(SNRc,[],4),[],5),1,1,1,n.Nexc,n.pH); 
    %max BiC/CO2 over Nexc, but min over pH values, take ratio (use for masking!)
for mm = 1:n.pH
    pHerr(:,:,:,:,mm) = abs(pHcalc(:,:,:,:,mm) - pH(mm)); %pH error  
    SNRc90(mm) = 1/(BtoC(mm)+1)/noise0; %SNR if CO2 was hit with one 90ï¿½ pulse
end

% Set variables for initial mask specification
%
pHtol = 0.05; %tolerable ï¿½error for pH accuracy
SNRbtocmin = 1; %constrains to SNRb/SNRc >= specified value when finding maxima
% SNRbtocmin = 1/5; %constrains to SNRb/SNRc >= specified value when finding maxima

SNRbmaskd = zeros(size(SNRb));
SNRcmaskd = SNRbmaskd;
SNRcmax = zeros(size(pH));
pHflg = true;%false; %if true, enforces pH accuracy constraint on data (can
    %set dynamically using figure buttons)
snrflg = true;%false; %if true, enforces SNR ratio constraint on data (can
    %set dynamically using figure buttons)

% Make pH accuracy + SNR mask, find max SNRc
%
makeMask; findMax;


%% INTERACTIVE FIGURE
%
% Prepare variables for interactive figure
%
parnames = [{'FAb'},{'FAc'},{'TR'},{'Nexc'}];
actXname = parnames{1}; %name of active x-variable on plot
actYname = parnames{2}; %name of active y-variable on plot
actpH = 1; %index of active pH value for plot
onames = parnames(~strcmp(parnames,actXname) & ~strcmp(parnames,actYname));
    %names of parameters currently not plotted
onamesX = parnames(~strcmp(parnames,actYname)); %names not chosen as x-axis
onamesY = parnames(~strcmp(parnames,actXname)); %names not chosen as x-axis
oval = ind.(onames{1})(1,actpH);
oval2 = ind.(onames{2})(1,actpH);
scrsz = get(groot,'ScreenSize');
sppi = get(groot,'ScreenPixelsPerInch');

% Create interactive figure for choosing parameter space to plot,
% restraints on pH accuracy + SNR
%
figure('Position',[0,scrsz(4)*1/6,scrsz(3)*2/3,scrsz(4)*5/6]*sppi,'Units','inches');
plotFig;
bg1 = uibuttongroup('Position',[0 0 1 .15]);
%axes selection
uicontrol(bg1,'Style','text','Position',[0 60 100 20],...
    'String','X Axis Parameter:');
xa = uicontrol(bg1,'Style','popupmenu','Position',[0 35 100 30],...
    'String',onamesX,'Callback',@xSelect);
uicontrol(bg1,'Style','text','Position',[0 20 100 20],...
    'String','Y Axis Parameter:');
ya = uicontrol(bg1,'Style','popupmenu','Position',[0 -5 100 30],...
    'String',onamesY,'Callback',@ySelect);
%parameter value selection
o1 = uicontrol(bg1,'Style','text','Position',[110 60 140 20],...
    'String',['Value of ' onames{1} ': ' num2str(pars.(onames{1})(oval))]);
s1 = uicontrol(bg1,'Style','slider','Min',1,'Max',n.(onames{1}),...
    'Value',ind.(onames{1})(1,actpH),'SliderStep',[1/n.(onames{1}) 3/n.(onames{1})],...
    'Position',[110 35 140 30],'Callback',@oValSet);
o2 = uicontrol(bg1,'Style','text','Position',[110 20 140 20],...
    'String',['Value of ' onames{2} ': ' num2str(pars.(onames{2})(oval2))]);
s2 = uicontrol(bg1,'Style','slider','Min',1,'Max',n.(onames{2}),...
    'Value',ind.(onames{2})(1,actpH),'SliderStep',[1/n.(onames{2}) 3/n.(onames{2})],...
    'Position',[110 -5 140 30],'Callback',@oVal2Set);
%set unplotted parameters to max values
uicontrol(bg1,'Style','pushbutton','Position',[260 40 120 30],...
    'String','<-Set for max SNR','Callback',@oValOpt);
uicontrol(bg1,'Style','pushbutton','Position',[260 0 120 30],...
    'String','<-Set for max SNR','Callback',@oVal2Opt);
%pH selection
uicontrol(bg1,'Style','text','Position',[400 45 80 20],...
    'String','pH:');
uicontrol(bg1,'Style','popupmenu','Position',[400 20 80 30],...
    'String',num2str(pH'),'Callback',@pHSelect);
%threshold constraints editing
uicontrol(bg1,'Style','text','Position',[520 45 80 15],...
    'String','Max |pH error|:');
uicontrol(bg1,'Style','edit','Position',[520 20 80 20],...
    'String',num2str(pHtol,'%1.2f'),'Callback',@setpHTol);
pb = uicontrol(bg1,'Style','pushbutton','Position',[520 0 100 20],...
    'String','Turn pH mask off','Callback',@selpHMask); 
uicontrol(bg1,'Style','text','Position',[640 45 80 15],...
    'String','Min SNRb/SNRc:');
uicontrol(bg1,'Style','edit','Position',[640 20 80 20],...
    'String',num2str(SNRbtocmin,'%2.3f'),'Callback',@setSNRRatio);
sb = uicontrol(bg1,'Style','pushbutton','Position',[640 0 100 20],...
    'String','Turn SNR mask off','Callback',@selSNRMask);


%% INTERNAL FUNCTIONS
%
% makeMask: Take values for pHtol and SNRbtocmin, mask SNR values based on
% these constraints
%
function makeMask
    mask = true(size(SNRc));
    if pHflg
        mask = mask .* (pHerr < pHtol); %pH accuracy constraint mask
        mask = repmat(prod(mask,5),[1 1 1 1 n.pH]); %enforces pH accuracy across all simulated pH values
    end
    if snrflg
        mask = mask .* (SNRbtocchk >= SNRbtocmin);
    end
    SNRbmaskd = SNRb .* mask;
    SNRcmaskd = SNRc .* mask;
end

% findMax: Finds the max SNRc for the specified constraints on pH and SNR,
% displays in main window
%
function findMax
    clear ind
    disp('Constraints:')  
    disp(['|pHerr|' tab 'minSNR'])
    disp([num2str(pHtol*pHflg,'%1.2f') tab num2str(SNRbtocmin*snrflg,'%1.3f')])
    for mm = 1:n.pH
        SNRcmax(mm) = max(reshape(SNRcmaskd(:,:,:,:,mm),1,[]));
        if SNRcmax(mm) ~= 0 && ~isnan(SNRcmax(mm)) %catch if SNRcmax is invalid value
            ctr = 0; %counts how many maxima found
            for ii = 1:n.FAb
                for jj = 1:n.FAc
                    for kk = 1:n.TR
                        for ll = 1:n.Nexc
                            if SNRcmaskd(ii,jj,kk,ll,mm) == SNRcmax(mm)
                                ctr = ctr + 1;
                                ind.FAb(ctr,mm) = ii;
                                ind.FAc(ctr,mm) = jj;
                                ind.TR(ctr,mm) = kk;
                                ind.Nexc(ctr,mm) = ll;
                            end
                        end
                    end
                end
            end
            disp(['pH = ' num2str(pH(mm),'%1.1f') ': ' num2str(size(ind.FAb,1),'%i') ...
                ' CO2 SNR maximum(a) found'])
            disp(['FAb(°)' tab 'FAc(°)' tab 'TR(ms)' tab 'Nexc' tab 'Ttot(s)' tab ...
                'SNRb/c' tab 'pHerr' tab 'Gain over 90°'])
            for ii = 1:size(ind.FAb,1)
                b = min(max(SNRbmaskd(ind.FAb(ii,mm),ind.FAc(ii,mm),ind.TR(ii,mm),:,:),[],4),[],5);
                    %use the minimum SNRb over all pH values, but maximized
                    %by Nexc
                c = SNRcmaskd(ind.FAb(ii,mm),ind.FAc(ii,mm),ind.TR(ii,mm),ind.Nexc(ii,mm),mm);
                disp([num2str(pars.FAb(ind.FAb(ii,mm)),'%2.0f') tab ...
                    num2str(pars.FAc(ind.FAc(ii,mm)),'%2.0f') tab ...
                    num2str(pars.TR(ind.TR(ii,mm)) * 1000,'%2.0f') tab ...
                    num2str(pars.Nexc(ind.Nexc(ii,mm)),'%2.0f') tab ...
                    num2str(pars.TR(ind.TR(ii,mm)) * pars.Nexc(ind.Nexc(ii,mm)),'%2.1f') tab ...
                    num2str(b/c,'%2.3f') tab ...
                    num2str(pHerr(ind.FAb(ii,mm),ind.FAc(ii,mm),ind.TR(ii,mm),ind.Nexc(ii,mm),mm),'%1.3f') tab ...
                    num2str(c/SNRc90(mm),'%2.3f')])
                % If maximum occurs for any min/max value of a parameter, display warning
                % message
                %
                if ind.FAb(ii,mm) == 1
                    warning('Max SNR occurred at min value for FAb. Consider re-running with lower min FAb')
                elseif ind.FAb(ii,mm) == n.FAb
                    warning('Max SNR occurred at max value for FAb. Consider re-running with higher max FAb')
                end
                if ind.FAc(ii,mm) == 1
                    warning('Max SNR occurred at min value for FAc. Consider re-running with lower min FAc')
                elseif ind.FAc(ii,mm) == n.FAc
                    warning('Max SNR occurred at max value for FAc. Consider re-running with higher max FAc')
                end
                if ind.TR(ii,mm) == 1
                    warning('Max SNR occurred at min value for TR. Consider re-running with lower min TR')
                elseif ind.TR(ii,mm) == n.TR
                    warning('Max SNR occurred at max value for TR. Consider re-running with higher max TR')
                end   
                if ind.Nexc(ii,mm) == 1
                    warning('Max SNR occurred at min value for Nexc. Consider re-running with lower min Nexc')
                elseif ind.Nexc(ii,mm) == n.Nexc
                    warning('Max SNR occurred at max value for Nexc. Consider re-running with higher max Nexc')
                end
            end
        end
    end
    disp('NOTE: SNRb/c relates MINIMUM SNRc to the MINIMUM SNRb over specified pH values, both maximized for Nexc')
end

% xSelect: Set parameter to plot along x-axis
%
function xSelect(source,~)
    actXname = onamesX{source.Value};
    onames = parnames(~strcmp(parnames,actXname) & ~strcmp(parnames,actYname));
    %reset values for other variables
    oval = ind.(onames{1})(1,actpH); 
    oval2 = ind.(onames{2})(1,actpH);
    updateGUI;
    plotFig;
end

% ySelect: Set parameter to plot along y-axis
%
function ySelect(source,~)
    actYname = onamesY{source.Value};
    onames = parnames(~strcmp(parnames,actXname) & ~strcmp(parnames,actYname));
    %reset values for other variables
    oval = ind.(onames{1})(1,actpH);
    oval2 = ind.(onames{2})(1,actpH);
    updateGUI;
    plotFig;
end

% pHSelect: Set pH value for plotting
%
function pHSelect(source,~)
    actpH = round(source.Value);
    plotFig;
end

% oValSet: Sets value of 1st unplotted variable, displays
%
function oValSet(source,~)
    oval = round(source.Value);
    set(o1,'String',['Value of ' onames{1} ': ' num2str(pars.(onames{1})(oval))]);
    plotFig;
end

% oVal2Set: Sets value of 2nd unplotted variable, displays
%
function oVal2Set(source,~)
    oval2 = round(source.Value);
    set(o2,'String',['Value of ' onames{2} ': ' num2str(pars.(onames{2})(oval2))]);
    plotFig;
end

% oValOpt: Sets value of 1st unplotted variable for max SNR
%
function oValOpt(~,~)
    oval = ind.(onames{1})(1,actpH);
    set(o1,'String',['Value of ' onames{1} ': ' num2str(pars.(onames{1})(oval))]);
    set(s1,'Value',oval);
    plotFig;
end

% oVal2Opt: Sets value of 2nd unplotted variable for max SNR
%
function oVal2Opt(~,~)
    oval2 = ind.(onames{2})(1,actpH);
    set(o2,'String',['Value of ' onames{2} ': ' num2str(pars.(onames{2})(oval2))]);
    set(s2,'Value',oval2);
    plotFig;
end

% setpHTol: sets value for maximum allowed |pH error|
%
function setpHTol(source,~)
    pHtol = str2double(source.String);
    makeMask; findMax;
    plotFig;
end

% selpHMask: Toggle pH accuracy masking on/off
%
function selpHMask(~,~)
    pHflg = ~pHflg;
    if pHflg %update button string based on value
        set(pb,'String','Turn pH mask off');
    else
        set(pb,'String','Turn pH mask on');
    end
    makeMask; findMax;
    plotFig;
end

% setSNRRatio: sets value for minimum allowed SNRb/SNRc
%
function setSNRRatio(source,~)
    SNRbtocmin = str2double(source.String);
    makeMask; findMax;
    plotFig;
end

% selSNRMask: Toggle pH accuracy masking on/off
%
function selSNRMask(~,~)
    snrflg = ~snrflg;
    if snrflg %update button string based on value
        set(sb,'String','Turn SNR mask off');
    else
        set(sb,'String','Turn SNR mask on');
    end
    makeMask; findMax;
    plotFig;
end

% updateGUI: Updates all values in interactive area of figure
%
function updateGUI
    set(o1,'String',['Value of ' onames{1} ': ' num2str(pars.(onames{1})(oval))]);
    set(o2,'String',['Value of ' onames{2} ': ' num2str(pars.(onames{2})(oval2))]);
    set(s1,'Min',1,'Max',n.(onames{1}),'Value',ind.(onames{1})(1,actpH),...
        'SliderStep',[1/n.(onames{1}) 3/n.(onames{1})]);
    set(s2,'Min',1,'Max',n.(onames{2}),'Value',ind.(onames{2})(1,actpH),...
        'SliderStep',[1/n.(onames{2}) 3/n.(onames{2})]);
    onamesX = parnames(~strcmp(parnames,actYname)); %names not chosen as y-axis
    onamesY = parnames(~strcmp(parnames,actXname)); %names not chosen as x-axis
    set(xa,'String',onamesX,'Value',find(strcmp(onamesX,actXname))); 
        %update selection options so you can't choose same parameter for x and y
    set(ya,'String',onamesY,'Value',find(strcmp(onamesY,actYname))); 
        %update selection options so you can't choose same parameter for x and y
end

% plotFig: Plots surface plot of SNRc vs the two selected parameters, at
% the parameter values specified for the other two parameters and at the 
% specified pH.
%
function plotFig
    %First, permute SNRcmaskd so that the dimensions are in the order
    %[Xpar, Ypar, 1stunusedpar, 2ndunusedpar, pH]
    indices(1) = find(strcmp(parnames,actYname));
    indices(2) = find(strcmp(parnames,actXname));
    indices(3) = find(strcmp(parnames,onames{1}));
    indices(4) = find(strcmp(parnames,onames{2}));
    plotSNR = permute(SNRcmaskd,[indices,5]);
    surf(pars.(actXname),pars.(actYname),squeeze(plotSNR(:,:,oval,oval2,actpH)))
    set(gca,'OuterPosition',[0 .2 1 .8]);
    axis([0 inf 0 inf 0 SNRcmax(actpH)]);
    axis('square');
    title(['pH ' num2str(pH(actpH),'%1.1f')])
    xlabel(actXname)
    ylabel(actYname)
    zlabel('CO_2 SNR, masked')    
end


end


%% EXTERNAL FUNCTIONS
%
% bmsim: Simulates z-magnetization exchange + T1 relaxation between pools A 
% and B over a given TR, using first-order forward and reverse exchange
% rates
%
function [Maf,Mbf] = bmsim(Ma0,Mb0,kab,kba,TR,T1,dt)
np = floor(TR / dt); %# of points for simulation
M = zeros (2 , np);
M(:,1) = [Ma0; Mb0];
R1a = 1 / T1;
R1b = R1a;
dtp = TR / (np-1);
% Matob = 0; %magnetization that goes from metabolite A to metabolite B
for i = 1 : np-1
     K = [-1*(R1a+kab)  kba             ;...
          kab           -1*(R1b+kba)    ];
     M(:,i+1) = expm(K * dtp) * M(:,i); %exponential model
%      Matob = Matob + M(1,i) * (exp(kab * dtp) - 1) - M(2,i) * (exp(kba * dtp) - 1); 
%         %net magnetization that went from A -> B over dtp
end
Maf = M(1,end);
Mbf = M(2,end);
end