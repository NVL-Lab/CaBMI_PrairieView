function [simulations_went_wrong] = run_posthoc_simulation(folder_list, experiment_types, number_ts)
%{
Function to run the simulations both ways from a table df that contains
info about the experiments (normally obtained posthoc with python analysis)
inputs examples:
folder_list = struct('FA', 'D:/data', 'FB', 'F:/data', 'FC', 'G:/data');
folder_df = 'C:/Users/Nuria/Documents/DATA/D1exp/df_data';
calibration_frames = 27000;
number_ts = 100;
experiment_types = {'D1act', 'CONTROL', 'CONTROL_LIGHT', 'CONTROL_AGO', 'RANDOM', 'NO_AUDIO', 'DELAY'}
%}

    %% initiate some vars
    [tset] = define_BMI_task_settings();
    [fbset]   = define_fb_audio_settings();
    frames_per_reward_range = tset.cb.sec_per_reward_range * tset.im.frameRate;
    simulations_went_wrong = {};
    
    %% iterate through all the experiments
    for e = 1:length(experiment_types)
        name = ['df_',  experiment_types{e},  '.parquet'];
        df = parquetread(fullfile(folder_df, name));
    
    %% iterate through all the mice
        for i = 1:size(df, 1)
            clear neurons;
            row = df(i, :);
            disp('Processing row: ' + row.session_path);
            % obtain the correct folders where the data is stored
            data_path = find_data_path(folder_list, row.mice_name);
            folder_process = fullfile(data_path, 'process', row.session_path);
            folder_suite2p = fullfile(folder_process, 'suite2p/plane0');
            folder_save = fullfile(folder_process, 'simulation_posthoc');
            if not(isfolder(folder_save))
                mkdir(folder_save)
            end
            dff = readNPY(fullfile(folder_suite2p, 'dff.npy'));
            is_cell = readNPY(fullfile(folder_suite2p, 'iscell.npy'));
            direct_neurons = load(fullfile(folder_suite2p, 'direct_neurons.mat'));
            indirect_neurons = find(is_cell(:,1)==1)-1;
            indirect_neurons = setdiff(indirect_neurons, [direct_neurons.E1,  direct_neurons.E2, ...
                direct_neurons.exclude]);
            neurons.T1 = direct_neurons;
            neurons.T2.E1 = direct_neurons.E2;
            neurons.T2.E2 = direct_neurons.E1;
            if length(indirect_neurons) >= 4
                for t=3:number_ts+2
                    randomIndices = randperm(length(indirect_neurons), 4);
                    neurons.(['T' num2str(t)]).E1 = indirect_neurons(randomIndices(1:2))';
                    neurons.(['T' num2str(t)]).E2 = indirect_neurons(randomIndices(3:4))';
                end
            end
            Ts = fieldnames(neurons);
            % obtain target T1
            
            for t = 1:numel(Ts)
                tname = Ts{t};
                tval = neurons.(tname);
                ensemble_neurons = [tval.E1  tval.E2];
                [target_info_path, ~, ~] = baseline2targetposthoc(dff(ensemble_neurons+1, 1:calibration_frames)', ...
                tval.E1 + 1, tval.E2 + 1, frames_per_reward_range, tset, folder_save, fbset, tname);
                close all
                target_info = load(target_info_path);

                % run simulation of T1 using same target_file     
                simulated_data = BMI_simulation_posthoc(dff(ensemble_neurons+1, calibration_frames:end), tset, target_info);
                close all

                data = simulated_data;
                bData = target_info;
                save(fullfile(folder_save, ['simulated_data_', tname, '.mat']), 'data', 'bData')
            end
        end
    end
