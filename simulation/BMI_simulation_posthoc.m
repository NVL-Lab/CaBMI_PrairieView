function [data] = BMI_simulation_posthoc(dff, tset, target_info)
%{
Function to simulate the BMI  
%}    
    
    %% BMI parameters 
    relaxationTime = 0;  % there can't be another hit in this many sec
    nonBufferUpdateCounter = tset.prefix_win;  %counter when we dont want to update the buffer: 
    
    % values of parameters in frames
    expectedLengthExperiment = length(dff(1,:)); % in frames
    relaxationFrames = round(relaxationTime * tset.im.frameRate);
    
    %% Define BMI parameters from baseline calibration
    back2Base = 1/2*target_info.T1; % cursor must be under this value to be able to hit again

    %% INITIALIZE
    %pre-allocating arrays
    Fbuffer = single(nan(length(target_info.E_id), tset.dff_win));  %define a windows buffer
    data.cursor = single(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.trialStart = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector trialStart
    data.selfHits = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits
    
    %initializing counters 
    data.selfTargetCounter = 0; 
    data.trialCounter = 0; %todo remove one
    BufferUpdateCounter = 0;
    back2BaseCounter = 0;
    
    %initializing flags
    trialFlag = 1;
    backtobaselineFlag = 0;
    
    %% SIMULATION 
    
    data.frame = 1;

    for i=1:expectedLengthExperiment
        if nonBufferUpdateCounter == 0
            Fbuffer(:, 1:end-1) = Fbuffer(:, 2:end);
            Fbuffer(:,end) = dff(:,data.frame);
            
            [~, cursor_i, target_hit, ~, ~, ~, ~, ~] = ...
                dff2cursor_target(single(nanmean(Fbuffer, 2)), target_info, tset.cursor_zscore_bool);
            data.cursor(data.frame) = cursor_i;

            if (BufferUpdateCounter == 0) 
                % Is it a new trial?
                if trialFlag && ~backtobaselineFlag
                    data.trialStart(data.frame) = 1;
                    data.trialCounter = data.trialCounter + 1;
                    trialFlag = 0;
                end

                if backtobaselineFlag 
                    if data.cursor(data.frame) <= back2Base 
                        back2BaseCounter = back2BaseCounter+1;
                    end
                    if back2BaseCounter >= tset.back2BaseFrameThresh
                        backtobaselineFlag = 0;
                        back2BaseCounter = 0;
                    end
                else
                    if target_hit      %if it hit the target
                        %Self hit:
                        data.selfTargetCounter = data.selfTargetCounter + 1;
                        data.selfHits(data.frame) = 1;

                        BufferUpdateCounter = relaxationFrames; 
                        backtobaselineFlag = 1;
                        % update trials and hits vector
                        trialFlag = 1;                                     
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

        data.frame = data.frame + 1;

    end   
end

