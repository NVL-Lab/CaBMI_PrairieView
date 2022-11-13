function [mat_path] = Baseline_Acqnvs_Prairie(path_data, roi_mask, tset, a, apin)
 %{
Function to acquire the baseline in a prairie scope
animal -> animal for the experiment
day -> day for the experiment
neuronMask -> matrix for spatial filters with px*py*unit 

%}

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    dilation_factor = 2; 
    expectedLengthExperiment = ...
        ceil(tset.cb.baseline_len*tset.im.frameRate*dilation_factor); 
    %in frames
    
    %save directory: 
    %mat file: 
    mat_path = ...
        fullfile(path_data.savePath, ['BaselineOnline' datestr(datetime('now'), 'yymmddTHHMMSS') '.mat']);
    

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************

    global pl baseActivity
    
    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(mat_path));  %in case of ctrl-c it will launch cleanmeup

    %% Prepare the nidaq
    s = daq.createSession('ni');
    addDigitalChannel(s,'dev6','Port0/Line0:2','OutputOnly');
    ni_out = [0 0 0]; 
    outputSingleScan(s,ni_out);%set   
    ni_getimage = [0 1 0]; 
    
    %% prepare the arduino
    % give random reward to trigger the jetball
    a.writeDigitalPin("D9", 1); pause(1);a.writeDigitalPin("D9",0)

    %% Prepare for Prairie
    % connection to Prairie
    pl = actxserver('PrairieLink.Application');
    pl.Connect()
    
    % pause needed for prairie to respond
    pause(2)

    % Prairie variables
    px = pl.PixelsPerLine();
    py = pl.LinesPerFrame();

    % Prairie commands
    pl.SendScriptCommands("-srd True 0");
    pl.SendScriptCommands("-lbs True 0");
    
    %define where to save the file
    savePathPrairie = fullfile(path_data.savePath, 'im');
    if ~exist(savePathPrairie, 'dir')
        mkdir(savePathPrairie);
    end
    savePrairieFiles(path_data.savePath, pl, 'baseline')

    lastFrame = zeros(px, py); % to compare with new incoming frames

    % set the environment for the Time Series in PrairieView
    loadCommand = "-tsl " + path_data.baseline_env;
    pl.SendScriptCommands(loadCommand);  
    
    %% Load and initialize Baseline variables

    % create smaller versions of the spatial filter

    numberNeurons = max(max(roi_mask));
    strcMask = obtain_Strc_Mask_from_Mask(roi_mask);    
    
    baseActivity = zeros(numberNeurons, expectedLengthExperiment) + nan;
    %% 
    %************************************************************************
    %*************************** RUN ********************************
    %************************************************************************
    frame = 1; % initialize frames
    %start the time_series scan
    %May need to put a break point on the next line, sometimes prairie
    %won't start scanning on it:
    pause(2); 
    pl.SendScriptCommands("-ts");  
    disp('sent -ts, pausing'); 
    pause(5);  %empirically discovered time for the prairie to start gears
    counterSame = 0;
    disp('Starting baseline acquisition')
    while counterSame < 500
        Im = pl.GetImage_2(tset.im.chan_data.chan_idx, px, py);
        if ~isequal(Im,lastFrame)   
            tic
            lastFrame = Im;   % comparison and assignment takes ~4ms
            outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,[0 0 0]);
            unitVals = obtain_Roi(Im, strcMask); % function to obtain Rois values 
            baseActivity(:,frame) = unitVals;
            frame = frame + 1;
            counterSame = 0;
            t = toc;
            if t < 1/(tset.im.frameRate*1.2)
                pause(1/(tset.im.frameRate*1.2) -t)
            end
        else
            counterSame = counterSame + 1;
        end

    end
    disp('Finished baseline')
    playTone(a, apin, 7000, 1)
   % pl.Disconnect();

end

function cleanMeUp(mat_path)
    global pl baseActivity
    disp('cleaning')
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(mat_path, 'baseActivity')
    if pl.Connected()
        pl.Disconnect();
    end
end





