global isOctave;
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if ~isdeployed
addpath(genpath(fileparts(mfilename('fullpath'))))


% Remove temp dir from path if it exists
tmpDir = fullfile(pwd, 'tmp');
if exist(tmpDir, 'dir')
    rmpath(genpath(tmpDir))
end


if ~isOctave % MATLAB
    % Check if specific toolboxes are installed
    if ~license('test', 'Optimization_Toolbox'), warning('Optimization Toolbox is not installed on your system: kinetic modeling won''t work. Please consider installing <a href="matlab:matlab.internal.language.introspective.showAddon(''OP'');">Optimization Toolbox</a>'); end
    if ~license('test', 'Image_Toolbox'), warning('Image Toolbox is not installed: some display functions will not work. Consider installing <a href="matlab:matlab.internal.language.introspective.showAddon(''IP'');">Image Processing Toolbox</a>'); end
    if ~license('test', 'Statistics_Toolbox'), warning('Statistics and Machine Learning Toolbox is not installed: some functions will not work.'); end
else % OCTAVE
    % install octave package
    installlist = {'optim','image','io','statistics', 'signal'};
    for ii=1:length(installlist)
        try
            disp(['loading ' installlist{ii}])
            pkg('load',installlist{ii})
        catch
            errorcount = 1;
            while errorcount % try to install 30 times (Travis)
                try
                    pkg('install','-forge',installlist{ii})
                    pkg('load',installlist{ii})
                    errorcount = 0;
                catch err
                    errorcount = errorcount+1;
                    if errorcount>30
                        error(err.message)
                    end
                end
            end
        end
    end

end
end