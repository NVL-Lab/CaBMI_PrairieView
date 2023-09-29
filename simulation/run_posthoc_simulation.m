function [simulations_went_wrong] = run_posthoc_simulation(folder_list, experiment_types, number_ts)
%{
Function to run the simulations both ways from a table df that contains
info about the experiments (normally obtained posthoc with python analysis)
inputs examples:
folder_list = struct('FA', 'D:/data', 'FB', 'F:/data', 'FC', 'G:/data');
folder_df = 'C:/Users/Nuria/Documents/DATA/D1exp/df_data';
calibration_frames = 27000;
number_ts = 100;
min_T = 0.7
experiment_types = {'D1act', 'CONTROL', 'CONTROL_LIGHT', 'CONTROL_AGO', 'RANDOM', 'NO_AUDIO', 'DELAY'}
%}

    %% initiate some vars
    [tset] = define_BMI_task_settings();
    frames_per_reward_range = tset.cb.sec_per_reward_range * tset.im.frameRate;
    
    %% iterate through all the experiments
    for e = 1:length(experiment_types)
        name = ['df_',  experiment_types{e},  '.parquet'];
        df = parquetread(fullfile(folder_df, name));
    
    %% iterate through all the mice
        for i = 1:size(df, 1)
            clear neurons;
            row = df(i, :);
            % obtain the correct folders where the data is stored
            data_path = find_data_path(folder_list, row.mice_name);
            folder_process = fullfile(data_path, 'process', row.session_path);
            folder_suite2p = fullfile(folder_process, 'suite2p/plane0');
            folder_save = fullfile(folder_process, 'simulation_posthoc');
            if not(isfolder(folder_save))
                mkdir(folder_save)
            end
            f_raw = readNPY(fullfile(folder_suite2p, 'F.npy'));
            bad_frames_file = fullfile(folder_suite2p, 'bad_frames.mat');
            if exist(bad_frames_file, 'file') > 0
                bad_frames = load(bad_frames_file);
                f_good = f_raw(:, ~bad_frames.bad_frames_bool);
                calibration_frames = calibration_frames - sum(bad_frames.bad_frames_bool(1:calibration_frames));
            else
                f_good = f_raw;
            end
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
                    neurons.(['T' num2str(t)]) = struct();
                end
            end
            Ts = fieldnames(neurons);
            iterat_Ts= 100;
            t = 1;
            while t <= numel(Ts) && iterat_Ts > 0
                disp('Processing row: ' + row.session_path + ' ' + Ts{t});
                tname = Ts{t};
                iterat = 500;
                if t>2
                    Tnotfound = true;
                    while Tnotfound && iterat > 0
                        iterat = iterat - 1;
                        randomIndices = randperm(length(indirect_neurons), 4);
                        neurons.(tname).E1 = indirect_neurons(randomIndices(1:2))';
                        neurons.(tname).E2 = indirect_neurons(randomIndices(3:4))';
                        tval = neurons.(tname);
                        [target_info_path, max_iter_achieved] = baseline2targetposthoc(f_good(:, 1:calibration_frames), ...
                        tval.E1 + 1, tval.E2 + 1, frames_per_reward_range, tset, folder_save, tname);
                        close all
                        if ~max_iter_achieved
                            target_info = load(target_info_path);
                            disp(['T: ', num2str(target_info.T1), ' iter: ', num2str(iterat)]);
                        else
                            continue
                        end
                        if target_info.T1 > min_T
                            Tnotfound = false;
                        end
                    end
                else
                    tval = neurons.(tname);
                    [target_info_path, max_iter_achieved] = baseline2targetposthoc(f_good(:, 1:calibration_frames), ...
                    tval.E1 + 1, tval.E2 + 1, frames_per_reward_range, tset, folder_save, tname);
                    target_info = load(target_info_path);
                end
                ensemble_neurons = [tval.E1  tval.E2];
                % run simulation of T1 using same target_file     
                simulated_data = BMI_simulation(f_good(ensemble_neurons+1, calibration_frames:end), tset, target_info);
                if t>2 && ~check_arr_stability(simulated_data.cursor)
                    disp('Cursor was not stable this T will not be saved')
                    delete(fullfile(folder_save, ['target_calibration_' tname '.mat']));
                    delete(fullfile(folder_save, ['BMI_target_info_' tname '.mat']));
                    iterat_Ts = iterat_Ts - 1;
                else
                    t = t + 1;
                    simulated_data.good_sim = ~max_iter_achieved && target_info.T1 > min_T;
                    data = simulated_data;
                    bData = target_info;
                    save(fullfile(folder_save, ['simulated_data_', tname, '.mat']), 'data', 'bData')
                end
                if iterat <=0
                    disp('Max iterations to find good indirect neurons pairs')
                    delete(fullfile(folder_save, ['target_calibration_' tname '.mat']));
                    delete(fullfile(folder_save, ['BMI_target_info_' tname '.mat']));
                    delete(fullfile(folder_save, ['simulated_data_' tname '.mat']));
                    break
                end
            end
        end
    end
