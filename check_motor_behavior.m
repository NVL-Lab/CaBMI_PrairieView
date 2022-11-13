function check_motor_behavior(a, path_data, tset, expt_str)
%{
    Function to check the behavior of the animal while also recording activity
%}

    global pl
    
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
    loadCommand = ['-tsl ' path_data.baseline_env];
    pl.SendScriptCommands(loadCommand);   

    % set the path where to store the imaging data -SetSavePath (-p) "path" ["addDateTime"]
    savePrairieFiles(path_data.savePath, pl, expt_str)  
    
    %% prepare the arduino
    % give random reward to trigger the jetball
    a.writeDigitalPin("D9", 1); pause(1);a.writeDigitalPin("D9",0)
    
    %% 
    %************************************************************************
    %*************************** RUN ********************************
    %************************************************************************
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
            lastFrame = Im;   % comparison and assignment takes ~4ms
            counterSame = 0;
        else
            counterSame = counterSame + 1;
        end

    end
    disp('Finished baseline')
    playTone(a, apin, 7000, 1)
   % pl.Disconnect();

    
    




