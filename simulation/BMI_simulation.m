function [data] = BMI_simulation(bmiAct, tset, target_info)
%{
Function to simulate the BMI  
%}    
    
    %% BMI parameters 
    relaxationTime = 0;  % there can't be another hit in this many sec
    nonBufferUpdateCounter = tset.prefix_win;  %counter when we dont want to update the buffer: 
    initFrameBase = nonBufferUpdateCounter + 1;
    
    % values of parameters in frames
    expectedLengthExperiment = length(bmiAct(1,:)); % in frames
    relaxationFrames = round(relaxationTime * tset.im.frameRate);
    
    %% Define BMI parameters from baseline calibration
    back2Base = 1/2*target_info.T1; % cursor must be under this value to be able to hit again
    numberNeurons = length(target_info.E_id);

    %% INITIALIZE
    %pre-allocating arrays
    Fbuffer = single(nan(numberNeurons, tset.dff_win));  %define a windows buffer
    
    data.cursor = single(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.baseVector = single(nan(numberNeurons,ceil(expectedLengthExperiment)));  %define a very long vector for cursor    
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
    baseBuffer_full = 0; %bool indicating the Fbuffer filled
    
    %% SIMULATION 
    
    data.frame = 1;
    

    for i=1:expectedLengthExperiment
        if nonBufferUpdateCounter == 0
            
            % obtain value of the neurons fluorescence
            unitVals = bmiAct(:,data.frame);
            
            % update F buffer
            Fbuffer(:, 1:end-1) = Fbuffer(:, 2:end);
            Fbuffer(:,end) = unitVals;

            % calculate F0 baseline activity 
            if data.frame == initFrameBase
                baseval = single(ones(numberNeurons,1)).*unitVals/tset.f0_win;
            elseif ~baseBuffer_full && data.frame <= (initFrameBase+tset.f0_win)
                baseval = baseval + unitVals/tset.f0_win;
                if data.frame == (initFrameBase+tset.f0_win)
                    baseBuffer_full = 1;
                end
            else 
                baseval = (baseval*(tset.f0_win - 1) + unitVals)./tset.f0_win;
            end
            data.baseVector(:,data.frame) = baseval;

            %Smooth F
            Fsmooth = single(nanmean(Fbuffer, 2));

            if baseBuffer_full
                % calculate (smoothed) DFF
                dff = (Fsmooth - baseval) ./ baseval;
                %Passing smoothed dff to "decoder"
                [~, cursor_i, target_hit, ~, ~, ~, ~, ~] = ...
                    dff2cursor_target(dff, target_info, tset.cursor_zscore_bool);
                data.cursor(data.frame) = cursor_i;
            end

            if (BufferUpdateCounter == 0) && baseBuffer_full
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
    
    data.bmiAct = bmiAct;
    
end

