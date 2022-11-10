function [flagBMI, flagDRstim, flagStimRandom, flagWater] = ...
    expt2bmi_flags(expt_str)
%{
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
    
    expt_cell = {...
        'BMI', ...
        'HoloVTA_pretrain', ...
        'Holo_pretrain', ...
        'VTA_pretrain'}; 
    %}

    if(strcmp(expt_str, 'BMI_stim'))
        flagBMI         = true; 
        flagDRstim      = true; 
        flagStimRandom  = false;
        flagWater       = false;    
    elseif(strcmp(expt_str, 'RandomDRstim'))
        flagBMI         = true; 
        flagDRstim      = false; 
        flagStimRandom  = true;
        flagWater       = false;   
    elseif(strcmp(expt_str, 'no_stim'))
        flagBMI         = true; 
        flagDRstim      = false; 
        flagStimRandom  = false;
        flagWater       = false;   
    elseif(strcmp(expt_str, 'BMI_no_stim_water'))
        flagBMI         = true; 
        flagDRstim      = false; 
        flagStimRandom  = false;
        flagWater       = true;    
    elseif(strcmp(expt_str, 'BMI_stim_water'))
        flagBMI         = true; 
        flagDRstim      = true; 
        flagStimRandom  = false;
        flagWater       = true;  
    elseif(strcmp(expt_str, 'BMI_random_stim_water'))
        flagBMI         = true; 
        flagDRstim      = false; 
        flagStimRandom  = true;
        flagWater       = true;
    else
        disp('You did not enter an acceptable expt str')
    end