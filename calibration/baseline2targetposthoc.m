function [target_info_path, max_iter_achieved] = baseline2targetposthoc(f_calib,   ...
    E1_base, E2_base, frames_per_reward_range, tset, savePath, tname)
    %{
    inputs:
    - f_calib - contains matrix, neural fluorescence from baseline file, num_samples X num_neurons_baseline 
    -E1_base - E1 idxs in baseline data
    -E2_base - E2 idxs in baseline data
    -frames_per_reward_range - a range on how many frames should elapse before a
    reward is expected.  Used to calibrate the target patterns.
    -prefix_win - number of frames at start of baseline that calibration should NEGLECT.
      (these starting frames are just ignored.)
    -f0_win_bool - if true, estimate f0 with a window of activity.  if false, estimate f0 using the full baseline, 
    -f0_win - number of frames to use for estimating f0.  (This code uses a
    rolling average, as used during the BMI)
    -dff_win_bool - whether to smooth dff in a window
    -dff_win - number of frames to use for smoothing dff
    -savePath - directory to save baseline results in in.
    -cursor_zscore_bool - if 1, neural activity is zscored before going into
    -cursor calculation. if 0, neural activity is not zscored.  
    -f0_init_slide - if 0, f0 is only used after f0_win samples.  if 1, f0 is
    adapted in the window from 0 to f0_win samples.
    -E2mE1_prctile: it is the lowest acceptable E2mE1_prctile for deciding the
    target threshold.  
    %}


    %1) Select Temporal: BMI E data from baseline 
    f_calib(:,isnan(f_calib(1,:))) = [];
    f_calib = f_calib.'; %num_samples X num_neurons

    %Throw out prefix frames:
    E1_raw = f_calib((tset.prefix_win+1):end,E1_base); 
    E2_raw = f_calib((tset.prefix_win+1):end,E2_base); 
    f_raw = [E1_raw E2_raw]; %first E1, then E2

    %2) Ensemble information
    num_E1 = length(E1_base); 
    num_E2 = length(E2_base); 
    num_neurons = num_E1 + num_E2;

    E_id = [1*ones(num_E1, 1); 2*ones(num_E2, 1)]; 
    E1_sel = E_id==1; 
    E1_sel_idxs = find(E1_sel); 
    E2_sel = E_id==2; 
    E2_sel_idxs = find(E2_sel); 



    %%
    %4) Decoder information
    %Decoder information

    [decoder, ~, ~, ~, ~] = def_decoder(num_neurons, E1_sel, E2_sel);

    %% obtain baseline of fluorescence f0
    %First process f0: 
    if(tset.cb.f0_win_bool)    
        %Calculate f0 as in BMI: 

        if tset.cb.f0_init_slide
            num_samples = size(f_raw,1);
            f0 = zeros(num_samples,num_neurons); 
            for i=1:length(f0)
                if i==1
                    f0(i,:) = f_raw(i,:);
                elseif i<tset.f0_win
                    f0(i,:) = (f0(i-1,:)*(i-1)+f_raw(i,:))/i;
                else
                    f0(i,:) = (f0(i-1,:)*(tset.f0_win-1)+f_raw(i,:))/tset.f0_win;
                end
            end
            f_postf0 = f_raw;
            f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);
        else
            num_samples = size(f_raw,1);
            f0 = zeros(num_samples-tset.f0_win+1, num_neurons); 
            f0(1,:) = mean(f_raw(1:tset.f0_win, :), 1);
            for i = 2:length(f0)
                f0(i,:) = f0(i-1,:)*((tset.f0_win-1)/tset.f0_win) + f_raw((i+tset.f0_win-1), :)/tset.f0_win; 
            end
            %Truncate data based on the tset.f0_win:
            f_postf0 = f_raw(tset.f0_win:end, :); 
            f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);        
        end
    else
        f_postf0 = f_raw; 
        f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);
        f0 = f0_mean; 
    end


    %% Compute dff and dff_z:
    
    %smooth f:
    if(tset.cb.dff_win_bool)
        num_samples = size(f_postf0,1);     
        f_smooth = zeros(num_samples, num_neurons); 
        smooth_filt = ones(tset.dff_win,1)/tset.dff_win;     
        for i=1:num_neurons
            f_smooth(:,i) = conv(f_postf0(:,i), smooth_filt, 'same'); 
        end
    else
        f_smooth = f_postf0; 
    end

    dff = (f_smooth-f0)./f0;

    %mean center the dff:
    n_mean = nanmean(dff,1); %1 x num_neurons
    mean_mat = repmat(n_mean, size(dff,1), 1);
    dffc = dff-mean_mat;
    %divide by std:
    n_std = nanstd(dffc, 0, 1); %var(dffc, 0, 1).^(1/2); %1 x num_neurons
    dff_z = dffc./repmat(n_std, [size(dff,1) 1]); 


    %%
    if tset.cursor_zscore_bool
        n_analyze = dff_z;
    else
        n_analyze = dff;
    end

    valid_idxs  = find(~isnan(n_analyze(:,1)));
    n_analyze   = n_analyze(valid_idxs, :); 
    analyze_cov = cov(n_analyze);
    analyze_mean = nanmean(n_analyze); %takes mean along dim 1.  n_analyze is num_samples X num_neurons.


    %%
    %Inputs: 
    %frames_per_reward_range
    %cov_bool
    reward_per_frame_range  = 1./frames_per_reward_range;
    E1_mean                 = mean(analyze_mean(E1_sel));
    E1_std                  = sqrt((E1_sel/num_E1)'*analyze_cov*(E1_sel/num_E1));
    E2_subord_mean          = zeros(num_E2,1);
    E2_subord_std           = zeros(num_E2,1); 
    for E2_i = 1:num_E2
        subord_sel                      = E2_sel;
        subord_sel(E2_sel_idxs(E2_i))   = 0; 
        E2_subord_mean(E2_i)            = mean(analyze_mean(subord_sel));     
        var_i                           = subord_sel'*analyze_cov*subord_sel; 
        E2_subord_std(E2_i)             = sqrt(var_i);     
    end

    %signals needed for target detection:
    cursor_obs                      = n_analyze*decoder; 

    %%
    %Iterate on T value, until perc correct value is achieved using truncated
    %neural activity

    %T:
    T0                  = max(cursor_obs);
    T                   = T0; 

    %E2:
    E2_coeff0           = 0.5;
    E2_coeff            = E2_coeff0; %multiplies the std dev, for figuring out E2_subord_thresh
    E2_subord_thresh    = E2_subord_mean+E2_subord_std*E2_coeff;

    %E1:
    E1_coeff0           = 0;
    E1_coeff            = E1_coeff0;
    E1_thresh           = E1_mean + E1_coeff*E1_std; %E1_mean_max; %E1_mean;


    T_delta = 0.05;
    task_complete = false;
    max_iter_achieved = false;
    max_iter = 10000;
    iter = 0;

%%

    while(~task_complete)

        %1) E2-E1 > alpha
        c1 = find(cursor_obs >= T); 
        %2) E1 < mu
        %c2 = find(E1_mean_analyze <= E1_thresh);
        %3) E2_subord > mu (anded with previous constraint)
        %For each idx, subtract the 
        %c3 = find(E2_subord_mean_analyze >= E2_subord_thresh(E2_dom_sel)); 
        hit_idxs_no_b2base = c1; %intersect(intersect(c1, c2), c3);
        %Remove hits that fall in a back2base

        %----------------------------------------------------------------------
        %Remove hits that fall in a back2base
        b2base_thresh = 0.5*T;
        hits_valid = ones(length(hit_idxs_no_b2base),1); 
        if length(hit_idxs_no_b2base) > 1
            for i = 2:length(hit_idxs_no_b2base)
                b2base_bool = sum(cursor_obs(hit_idxs_no_b2base(i-1):hit_idxs_no_b2base(i)) <= b2base_thresh) >= tset.back2BaseFrameThresh; 
                hits_valid(i) = b2base_bool; 
            end
        end
        hit_idxs_b2base         = hit_idxs_no_b2base(find(hits_valid)); 
        valid_hit_idxs          = hit_idxs_b2base;
        reward_prob_per_frame   = sum(hits_valid)/length(n_analyze);    
        %----------------------------------------------------------------------

        %Update T:
        if((reward_prob_per_frame >= reward_per_frame_range(1)) && (reward_prob_per_frame <= reward_per_frame_range(2)))
            task_complete = 1;
        elseif(reward_prob_per_frame > reward_per_frame_range(2))
            %Task too easy, make T harder:
            T = T+T_delta; 
        elseif(reward_prob_per_frame < reward_per_frame_range(1))
            %Task too hard, make T easier:
            T = T-T_delta; 

        end

        iter = iter+1;
        if(iter == max_iter)
            task_complete = true; 
            max_iter_achieved = true;
        end
    end



    %%
    num_hits_no_b2base = length(hit_idxs_no_b2base);
    num_valid_hits = length(valid_hit_idxs);
%
    %%
    %Save the results: 
    %1) All the steps here
    %2) Just the target parameters for running BMI

    %1)All the steps here
    clear h
    save_path = fullfile(savePath, ['target_calibration_' tname '.mat']); 
    save(save_path, 'num_hits_no_b2base', 'num_valid_hits'); 

    %2)Just the target parameters for running BMI
    target_info_file = ['BMI_target_info_' tname '.mat'];
    save_path = fullfile(savePath, target_info_file); 
    target_info_path = save_path; 
    %Change variable names for BMI code:
    T1 = T; %Change to T1, as this is what BMI expects
    save(save_path, 'n_mean', 'n_std', 'decoder', 'E_id', 'E1_sel_idxs', 'E2_sel_idxs', 'E1_base', 'E2_base', 'T1', 'E1_thresh', 'E1_coeff', 'E1_std', 'E2_subord_thresh', 'E2_coeff', 'E2_subord_mean', 'E2_subord_std'); 
end

function [decoder, E1_proj, E2_proj, E1_norm, E2_norm] = ...
    def_decoder(num_neurons, E1_sel, E2_sel)

    E1_proj = zeros(num_neurons, 1); 
    E1_proj(E1_sel) = 1;
    E1_norm = sum(E1_sel); %can replace with vector norm.  
    E1_proj = E1_proj/E1_norm;

    E2_proj = zeros(num_neurons, 1); 
    E2_proj(E2_sel) = 1; 
    E2_norm = sum(E2_sel); 
    E2_proj = E2_proj/E2_norm;

    decoder = E2_proj - E1_proj;
end
