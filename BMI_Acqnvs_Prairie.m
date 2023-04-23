function BMI_Acqnvs_Prairie(path_data, expt_str, baselineCalibrationFile, tset, vector_stim, ...
    debug_bool, debug_input, baseValSeed, fb_bool, fb_cal, a)
    %{
    Function to acquire the BMI in a prairie scope
    animal -> animal for the experiment
    day -> day for the experiment
    
    debug_input: num_neurons x num_samples

    baselineCalibrationFile:
    %AComp_BMI: matrix for spatial filters with px*py*unit
    %zscore parameters: n_mean, n_std
    %cursor parameters: decoder
    %target parameters: 'T1', 'E1_thresh', 'E2_subord_thresh'

    vector_stim -> vector of scheduled stims

    Flag description:
    flagBMI: addDigitalChannel
        determines if code detects self-generated target neural patterns
    flagDRstim: 
        determines if target neural patterns triggers DR stim
    flagStimRandom: 
        determines if DR stim will be delivered on a random schedule
    
    -fb_bool - bool to play tones
    -fb_cal - calibration for mapping cursor to audio feedback
    -a - arduino object for playing feedback tones.

    expt_str --> Experiments:
    0) BMI_stim
    flagBMI = true; (use self-generated hits to perform actions 
        for now actions are just to send DR stim
    flagDRstim = true; 
    flagStimRandom = false;
    flagWater = false;

    1) RandomDRstim
    flagDRstim = False; flagStimRandom = true;
    This follows the vector_stim schedule. 
    flagBMI = false; 
    flagWater = false

    2) no_stim
    flagBMI = true;
    flagDRstim = false; 
    flagStimRandom = false;
    flagWater = false;
    
    3) BMI_water
    flagBMI = true; (use self-generated hits to perform actions 
        for now actions are just to send DR stim
    flagDRstim = false; 
    flagStimRandom = false;
    flagWater = true;
    
    
    4) BMI_water_stim
    flagBMI = true; (use self-generated hits to perform actions 
        for now actions are just to send DR stim
    flagDRstim = true; 
    flagStimRandom = false;
    flagWater = true;
    
    5)BMI_water_random_stim
    flagBMI = true; (use self-generated hits to perform actions 
        for now actions are just to send DR stim
    flagDRstim = false; 
    flagStimRandom = true;
    flagWater = true;
    
%}

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    
    %% experiment FLAGS
    
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

    [flagBMI, flagDRstim, flagStimRandom, flagWater] = ...
        expt2bmi_flags(expt_str);
  
    
    %% BMI parameters 
    %tset.im.frameRate = 30; % TODO check if it can be obtained from prairie 
    relaxationTime = 0;  % there can't be another hit in this many sec
    %back2Base = 1/2*bData.T1; % cursor must be under this value to be able to hit again

    % values of parameters in frames
    expectedLengthExperiment = 60*30*tset.im.frameRate; % in frames
    relaxationFrames = round(relaxationTime * tset.im.frameRate);
    
    %% Load BMI parameters from baseline calibration
    bData = load(fullfile(baselineCalibrationFile));
    back2Base = 1/2*bData.T1; % cursor must be under this value to be able to hit again
%     single_bool = 1; 

    %Fields: 
    %'n_mean', 'n_std',
    %'AComp_BMI', 'T1', 'decoder', 'E_id', 'E1_sel_idxs', 'E2_sel_idxs', 
    %'E1_base', 'E2_base',
    %'E2_subord_thresh', 'E1_thresh', 'E2_coeff', 'E2_subord_mean', 'E2_subord_std'

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************
    
    global pl data
%     cursor hits trialStart bmiAct baseVector timeVector %TODO remove timeVector
    
    numberNeurons = length(bData.E_id);
    
    %pre-allocating arrays
    Fbuffer = single(nan(numberNeurons, tset.dff_win));  %define a windows buffer
    
    data.cursor = single(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.fb_freq    = single(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for fb_freq
    data.bmiAct = single(nan(numberNeurons, ceil(expectedLengthExperiment)));
%     data.bmidffz = single(nan(numberNeurons, ceil(expectedLengthExperiment)));
    data.baseVector = single(nan(numberNeurons,ceil(expectedLengthExperiment)));  %define a very long vector for cursor    
    data.selfHits = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits
    data.selfDRstim = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.vector_stim = vector_stim;
    data.vectorWater = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for water    
    data.randomDRstim = single(zeros(1,ceil(expectedLengthExperiment)));
    data.trialStart = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector trialStart
    %to debug!!! TODO REMOVE after debugging
    data.timeVector = single(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    
    if(debug_bool)
        data.fsmooth    = single(nan(numberNeurons, ceil(expectedLengthExperiment)));
        data.dff        = single(nan(numberNeurons, ceil(expectedLengthExperiment)));
        data.c1_bool    = single(nan(1, ceil(expectedLengthExperiment)));
        data.c2_val     = single(nan(1, ceil(expectedLengthExperiment)));
        data.c2_bool    = single(nan(1, ceil(expectedLengthExperiment)));
        data.c3_val     = single(nan(1, ceil(expectedLengthExperiment)));
        data.c3_bool    = single(nan(1, ceil(expectedLengthExperiment)));        
    end

    %initializing general flags and counters 
    data.selfTargetCounter = 0; 
    data.selfTarget_DR_stim_Counter = 0; 
    data.sched_random_stim = 0; 
    data.Water_Counter = 0;
     
    data.trialCounter = 0; %todo remove one
    trialFlag = 1;
    nonBufferUpdateCounter = tset.prefix_win;  %counter when we dont want to update the buffer: 
    initFrameBase = nonBufferUpdateCounter + 1;
    %beginning of experiment and VTA stim
    BufferUpdateCounter = 0;

    deliver_water       = 0;
    deliver_stim        = 0;
 
    
    back2BaseCounter = 0;
    backtobaselineFlag = 0;
    
    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(path_data.savePath, bData, debug_bool));  %in case of ctrl-c it will launch cleanmeup

%     %% Prepare the nidaq
    if(~debug_bool)
        clear s
        s = daq.createSession('ni');
        addDigitalChannel(s,'dev6','Port0/Line0:2','OutputOnly');
        ni_out = [0 0 0]; 
        outputSingleScan(s,ni_out);%set   
        ni_getimage = [0 1 0]; 
    end
%       Line 1: GetImage Pulse


    %% Prepare for Prairie
    % connection to Prairie
    if(~debug_bool)
        pl = actxserver('PrairieLink.Application');
        pl.Connect()
        pause(2);  % pause is needed to give time to Prairie to connect

        % Prairie variables
        px = pl.PixelsPerLine();
        py = pl.LinesPerFrame();

        % Prairie commands
        pl.SendScriptCommands('-srd True 0');
        pl.SendScriptCommands('-lbs True 0');

        % set the environment for the Time Series in PrairieView
        loadCommand = ['-tsl ' path_data.bmi_env];
        pl.SendScriptCommands(loadCommand);   

        % set the path where to store the imaging data -SetSavePath (-p) "path" ["addDateTime"]
        savePrairieFiles(path_data.savePath, pl, expt_str)  
    else
        px = 512; 
        py = 512; 
    end
    
    lastFrame = zeros(px, py); 
    %% load masks
    if(~debug_bool)
        load(fullfile(path_data.savePath, 'strcMask.mat'), 'strcMask');
    end

    %%
    
    %% ************************************************************************
    %*************************** RUN ********************************
    %************************************************************************

    %start the time_series scan
    if(~debug_bool)
        pause(2); 
        pl.SendScriptCommands('-ts');  
        pause(2);  %empirically discovered time for the prairie to start gears
    end
    data.frame = 1;
    
    % give random reward to trigger the jetball
    a.writeDigitalPin("D9", 1); pause(1);a.writeDigitalPin("D9",0)
    
    
    disp('STARTING RECORDING!!!')
    counterSame = 0; %Counts how many frames are the same as past,
    counterSameThresh = 500;
    baseBuffer_full = 0; %bool indicating the Fbuffer filled
    %---
    disp('baseBuffer filling!...')
    while (~debug_bool && counterSame < counterSameThresh) || (debug_bool && data.frame <= size(debug_input,2)) %while data.frame <= expectedLengthExperiment
        if ~debug_bool
            Im = pl.GetImage_2(tset.im.chan_data.chan_idx, px, py);
        else
            Im = zeros(px, py); 
        end
        
        if ~isequal(Im,lastFrame) || debug_bool
            tic; %Start timing to see length of an iteration
            if(~debug_bool)
                lastFrame = Im;   % comparison and assignment takes ~4ms
                % defines that we got an image
                outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,[0 0 0]);
            end
            
%             if nonBufferUpdat
            if nonBufferUpdateCounter == 0
                % obtain value of the neurons fluorescence
                if(~debug_bool)
                    unitVals = obtain_Roi(Im, strcMask); % function to obtain Rois values
                else
                    unitVals = debug_input(:,data.frame); 
                end
                data.bmiAct(:,data.frame) = unitVals;
                
                % update F buffer
                Fbuffer(:, 1:end-1) = Fbuffer(:, 2:end);
                Fbuffer(:,end) = unitVals;
                
                % calculate F0 baseline activity 
                if data.frame == initFrameBase
                    if ~isnan(sum(baseValSeed))
                        baseBuffer_full = 1; 
                        baseval = baseValSeed; 
                        disp('baseBuffer seeded!'); 
                    else
                        baseval = single(ones(numberNeurons,1)).*unitVals/tset.f0_win;
                    end
                    %---
                elseif ~baseBuffer_full && data.frame <= (initFrameBase+tset.f0_win)
%                     baseval = base(baseval*(data.frame - 1) + signal)./data.frame;
                    baseval = baseval + unitVals/tset.f0_win;
%                     disp(data.frame);
                    if data.frame == (initFrameBase+tset.f0_win)
                        baseBuffer_full = 1;
                        disp('baseBuffer FULL!'); 
                    end
                else %data.frame > (initFrameBase+tset.f0_win)
                    baseval = (baseval*(tset.f0_win - 1) + unitVals)./tset.f0_win;
                end
                data.baseVector(:,data.frame) = baseval;
                
                %Smooth F
                Fsmooth = single(nanmean(Fbuffer, 2));
                
                if debug_bool
                    data.fsmooth(:,data.frame) = Fsmooth;
                end
                
                if baseBuffer_full
                    %----------------------------------------------------------
                    %Cursor
                    % calculate (smoothed) DFF
                    dff = (Fsmooth - baseval) ./ baseval;
                    %Passing smoothed dff to "decoder"
                    [~, cursor_i, target_hit, c1_bool, c2_val, c2_bool, c3_val, c3_bool] = ...
                        dff2cursor_target(dff, bData, tset.cursor_zscore_bool);
%                     data.bmidffz(:,data.frame) = dff_z;
%--------------------------------------------------------------------------
                    disp(['Cursor: ' num2str(cursor_i)]); 
                    data.cursor(data.frame) = cursor_i;
%--------------------------------------------------------------------------                    
                    if debug_bool
                        data.dff(:,data.frame)      = dff;
                        data.c1_bool(data.frame)    = c1_bool; 
                        data.c2_val(data.frame)     = c2_val;
                        data.c2_bool(data.frame)    = c2_bool;
                        data.c3_val(data.frame)     = c3_val;
                        data.c3_bool(data.frame)    = c3_bool;
                    end
                    
                    %fb: 
%--------------------------------------------------------------------------                    
                    fb_freq_i = cursor2audio(cursor_i, fb_cal);  
%                     if(debug_bool)
%                         disp(['FB Freq: ' num2str(fb_freq_i)]);
%                     end
                    data.fb_freq(data.frame) = fb_freq_i;
%--------------------------------------------------------------------------

                    if(fb_bool && ~debug_bool)
                        %Send tone arduino
                        playTone(a,...
                            fb_cal.settings.arduino.pin,...
                            fb_freq_i,...
                            fb_cal.settings.arduino.duration)
                    end                    
                    
%                     disp(['Target : ' num2str(target_hit)]); 
%                     disp(['C1 - cursor: ' num2str(c1_bool)]); 
%                     disp(['C2 - E1 : ' num2str(c2_bool)]); 
%                     disp(['C3 - E2 subord : ' num2str(c3_bool)]);                     
                    % c1: cursor
                    % c2: E1_mean > E1_thresh
                    % c3: E2_subord_mean > E2_subord_thresh                    
                    %----------------------------------------------------------
                end
                
                if (BufferUpdateCounter == 0) && baseBuffer_full
%                     disp('HERE'); 
                    % Is it a new trial?
                    if trialFlag && ~backtobaselineFlag
                        data.trialStart(data.frame) = 1;
                        data.trialCounter = data.trialCounter + 1;
                        trialFlag = 0;
                        %start running the timer again
                        disp('New Trial!')
                    end

                    if backtobaselineFlag 
                        if data.cursor(data.frame) <= back2Base 
                            back2BaseCounter = back2BaseCounter+1;

                        end
                        if back2BaseCounter >= tset.back2BaseFrameThresh
                            backtobaselineFlag = 0;
                            back2BaseCounter = 0;
                            disp('back to baseline')
                        end
                    else
%                         disp('HERE2'); 
                        if target_hit      %if it hit the target
                            disp('target hit')
                            %Self hit:
                            data.selfTargetCounter = data.selfTargetCounter + 1;
                            data.selfHits(data.frame) = 1;
                            disp(['Trial: ', num2str(data.trialCounter), 'Num Self Hits: ', num2str(data.selfTargetCounter)]);

                            if(flagBMI)
                                
                                if(flagDRstim)
                                    deliver_stim = 1;
                                    data.selfTarget_DR_stim_Counter = data.selfTarget_DR_stim_Counter + 1;
                                    data.selfDRstim(data.frame) = 1;
                                    disp(['Trial: ', num2str(data.trialCounter), 'DR STIMS: ', num2str(data.selfTarget_DR_stim_Counter)]);  
                                end
                                if(flagWater)
                                    deliver_water = 1;
                                    data.Water_Counter = data.Water_Counter + 1;
                                    data.vectorWater(data.frame) = 1;
                                    disp(['Trial: ', num2str(data.trialCounter), 'Water: ', num2str(data.Water_Counter)]); 
                                end

                                disp('Target Achieved! (self-target)')
                                
                                if ~debug_bool
%                                         play(reward_sound);
                                    disp('RewardTone delivery!')
                                end                                        
                                BufferUpdateCounter = relaxationFrames; 
                                backtobaselineFlag = 1;
                                % update trials and hits vector
                                trialFlag = 1;                                     
                            end
                        end
                        if ~trialFlag
%                             disp(['HERE ' num2str(data.frame)]); 
                            if flagStimRandom
                                if ismember(data.frame, data.vector_stim)
                                    deliver_stim = 1;
                                    disp('SCHEDULED DR STIM'); 
                                    data.sched_random_stim = data.sched_random_stim + 1;
                                    %Also, save the frame we do this!!
                                    data.randomDRstim(data.frame) = 1;
                                end
                            end
                        end
                            
                    end
                else
                    if(BufferUpdateCounter>0)
                        BufferUpdateCounter = BufferUpdateCounter - 1;
                    end
                end
            else
                if(nonBufferUpdateCounter>0)
                    nonBufferUpdateCounter = nonBufferUpdateCounter - 1;
                end
            end
            
            if deliver_water 
                if(~debug_bool)
                    a.writeDigitalPin("D9", 1); pause(1);a.writeDigitalPin("D9",0)
                end
                deliver_water = 0; 
                disp('water delivered!'); 
            end
            if deliver_stim
                if tset.delay_flag
                    pause(tset.delay_time)
                end
                if(~debug_bool)
                    % blue light
                    a.writeDigitalPin("D5", 1); pause(0.2);a.writeDigitalPin("D5",0)
                    % UV light
                    a.writeDigitalPin("D3", 1); pause(1);a.writeDigitalPin("D3",0)
                end
                deliver_stim = 0; 
                disp('stim delivered!');
            end
                
            data.frame = data.frame + 1;
            data.timeVector(data.frame) = toc;
            counterSame = 0;
            if (~debug_bool && data.timeVector(data.frame) < 1/(tset.im.frameRate*1.2))
                pause(1/(tset.im.frameRate*1.2) - data.timeVector(data.frame))
            end
        else
            counterSame = counterSame + 1;
        end
    end
%    pl.Disconnect();
%     save(fullfile(savePath, ['BMI_online', datestr(datetime('now'), 'yymmddTHHMMSS'), '.mat']), 'data', 'bData')
end
% 
% % fires when main function terminates (normal, error or interruption)
function cleanMeUp(savePath, bData, debug_bool)
    global pl data
    disp('cleaning')
    % saving the global variables
    save(fullfile(savePath, ['BMI_online', datestr(datetime('now'), 'yymmddTHHMMSS'), '.mat']), 'data', 'bData')
    if ~debug_bool
        if pl.Connected()
            pl.Disconnect();
        end
    end
end

