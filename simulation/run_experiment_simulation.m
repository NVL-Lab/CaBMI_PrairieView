function [simulations_went_wrong] = run_experiment_simulation(df, folder_list)
%{
Function to run the simulations both ways from a table df that contains
info about the experiments (normally obtained posthoc with python analysis)
inputs examples:
folder_list = struct('FA', 'D:/data', 'FB', 'F:/data', 'FC', 'G:/data');
folder_df = 'C:/Users/Nuria/Documents/DATA/D1exp/df_data';
df_CONTROL = parquetread(fullfile(folder_df, 'df_CONTROL.parquet'));
df_CONTROL_AGO = parquetread(fullfile(folder_df, 'df_CONTROL_AGO.parquet'));
df_CONTROL_LIGHT = parquetread(fullfile(folder_df, 'df_CONTROL_LIGHT.parquet'));
df_D1act = parquetread(fullfile(folder_df, 'df_D1act.parquet'));
df_RANDOM = parquetread(fullfile(folder_df, 'df_RANDOM.parquet'));
df_DELAY = parquetread(fullfile(folder_df, 'df_DELAY.parquet'));
df_NO_AUDIO = parquetread(fullfile(folder_df, 'df_NO_AUDIO.parquet'));
df = df_D1act;
%}
    %% initiate some vars
    [tset] = define_BMI_task_settings();
    [fbset]   = define_fb_audio_settings();
    frames_per_reward_range = tset.cb.sec_per_reward_range * tset.im.frameRate;
    simulations_went_wrong = {};
    
    %% iterate through all the mice
    for i = 1:size(df, 1)
        row = df(i, :);
        disp('Processing row: ' + row.session_path);
        
        % obtain the correct folders where the data is stored
        data_path = find_data_path(folder_list, row.mice_name);
        folder_raw = fullfile(data_path, 'raw', row.session_path);
        folder_process = fullfile(data_path, 'process', row.session_path);
        
        n_f_file = fullfile(folder_raw, row.Baseline_online);
        roi_data_file = fullfile(folder_raw, row.roi_data);
        folder_save = fullfile(folder_process, 'simulation');
        if not(isfolder(folder_save))
            mkdir(folder_save)
        end
        
        % obtain and clean the raw data
        bmi_raw_data = clean_bmi_raw_data(fullfile(folder_raw, row.BMI_online));
        
        % run simulation of T1 using same target_file     
        simulated_data_T1 = BMI_simulation(bmi_raw_data.data.bmiAct, tset, bmi_raw_data.bData);
        if simulated_data_T1.selfTargetCounter ~= bmi_raw_data.data.selfTargetCounter
            disp('Something went wrong, simulated data is not the same as raw data');
            simulations_went_wrong = [simulations_went_wrong, row.session_path];
        end
        data = simulated_data_T1;
        bData = bmi_raw_data.bData;
        save(fullfile(folder_save, ['simulated_data_T1_', datestr(datetime('now'), 'yymmddTHHMMSS'), '.mat']), 'data', 'bData')
        % obtain a target_file for T2 based on baseline
        [target_info_path, ~, ~] = baseline2target(n_f_file, roi_data_file,  ...
            bmi_raw_data.bData.E2_base, bmi_raw_data.bData.E1_base, frames_per_reward_range, tset, folder_save, fbset);
        target_info_T2 = load(target_info_path);
        % run simulation of T2
        bmiAct_T2 = [bmi_raw_data.data.bmiAct(3:4, :);bmi_raw_data.data.bmiAct(1:2, :)];
        simulated_data_T2 = BMI_simulation(bmiAct_T2, tset, target_info_T2);
        data = simulated_data_T2;
        bData = target_info_T2;
        save(fullfile(folder_save, ['simulated_data_T2_', datestr(datetime('now'), 'yymmddTHHMMSS'), '.mat']), 'data', 'bData')
    end
