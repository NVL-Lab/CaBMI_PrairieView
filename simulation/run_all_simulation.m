%{
script to run all the simulations 

%}

folder_list = struct('FA', 'D:/data', 'FB', 'F:/data', 'FC', 'G:/data');
folder_df = 'C:/Users/Nuria/Documents/DATA/D1exp/df_data';
df_CONTROL = parquetread(fullfile(folder_df, 'df_CONTROL.parquet'));
df_CONTROL_AGO = parquetread(fullfile(folder_df, 'df_CONTROL_AGO.parquet'));
df_CONTROL_LIGHT = parquetread(fullfile(folder_df, 'df_CONTROL_LIGHT.parquet'));
df_D1act = parquetread(fullfile(folder_df, 'df_D1act.parquet'));
df_RANDOM = parquetread(fullfile(folder_df, 'df_RANDOM.parquet'));
df_DELAY = parquetread(fullfile(folder_df, 'df_DELAY.parquet'));
df_NO_AUDIO = parquetread(fullfile(folder_df, 'df_NO_AUDIO.parquet'));
simulations_wrong = {};

[s_wrong] = run_experiment_simulation(df_D1act, folder_list);
simulations_wrong = [simulations_wrong, s_wrong];
[s_wrong] = run_experiment_simulation(df_CONTROL, folder_list);
simulations_wrong = [simulations_wrong, s_wrong];
[s_wrong] = run_experiment_simulation(df_CONTROL_AGO, folder_list);
simulations_wrong = [simulations_wrong, s_wrong];
[s_wrong] = run_experiment_simulation(df_CONTROL_LIGHT, folder_list);
simulations_wrong = [simulations_wrong, s_wrong];
[s_wrong] = run_experiment_simulation(df_RANDOM, folder_list);
simulations_wrong = [simulations_wrong, s_wrong];
[s_wrong] = run_experiment_simulation(df_DELAY, folder_list);
simulations_wrong = [simulations_wrong, s_wrong];
[s_wrong] = run_experiment_simulation(df_NO_AUDIO, folder_list);
simulations_wrong = [simulations_wrong, s_wrong];

disp(simulations_wrong)


