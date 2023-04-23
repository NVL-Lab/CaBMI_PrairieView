function [a] = run_all_simulation(df)
%{
Function to run all the simulations both ways
inputs:
folder_process = 'F:\data\process'
folder_experiments = 'F:\data\raw'
df_BMI_CONTROL_AGO = parquetread(fullfile(folder_process, 'df_BMI_CONTROL_AGO.parquet'))
df_BMI_CONTROL_LIGHT = parquetread(fullfile(folder_process, 'df_BMI_CONTROL_LIGHT.parquet'))
df_BMI_CONTROL_RANDOM = parquetread(fullfile(folder_process, 'df_BMI_CONTROL_RANDOM.parquet'))
df_BMI_STIM_AGO = parquetread(fullfile(folder_process, 'df_BMI_STIM_AGO.parquet'))
df = df_BMI_STIM_AGO;
%}

    rows = height(df);
    [tset] = define_BMI_task_settings();
    frames_per_reward_range = tset.cb.sec_per_reward_range * tset.im.frameRate;
    for row = 1:rows 
        session_path = df(row, 'session_path').Variables;
        n_f_file = fullfile(folder_experiments, session_path, df(row, 'Baseline_online').Variables);
        roi_data_file = fullfile(folder_experiments, session_path, df(row, 'roi_data').Variables);
        strcMask_file = fullfile(folder_experiments, session_path, df(row, 'mask_data').Variables);
        load(fullfile(folder_experiments, session_path, df(row, 'BMI_online').Variables))
        folder_save = fullfile(folder_process, session_path, 'simulation');
        if not(isfolder(folder_save))
            mkdir(folder_save)
        end
        BMI_simulation(n_f_file, roi_data_file, strcMask_file, ...
            bData, frames_per_reward_range, tset, folder_save, data.cursor);
    end
