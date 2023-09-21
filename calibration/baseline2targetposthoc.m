function [target_info_path, target_cal_ALL_path, fb_cal] = baseline2targetposthoc(dff, ...
    E1_base, E2_base, frames_per_reward_range, tset, savePath, fbset, var_str)
%{
4.18.19
inputs:
-dff - contains matrix, neural fluorescence from baseline file, num_samples X num_neurons_baseline 
-roi_data_file - 
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

%%sa
%Select BMI ensemble temporal activity from baseline
%Create Ensemble Information
%Select BMI ensemble spatial components
%Decoder Information

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


%% prepare dff:
n_mean = nanmean(dff,1); %1 x num_neurons
mean_mat = repmat(n_mean, size(dff,1), 1);
dffc = dff-mean_mat;
%divide by std:
n_std = nanstd(dffc, 0, 1); %var(dffc, 0, 1).^(1/2); %1 x num_neurons

%%
 
valid_idxs  = find(~isnan(dff(:,1)));
dff   = dff(valid_idxs, :); 
analyze_cov = cov(dff);
analyze_mean = nanmean(dff); %takes mean along dim 1.  dff is num_samples X num_neurons.


%% obtain dff as in the bmi
Fbuffer = single(nan(size(dff,1), tset.dff_win));  %define a windows buffer
dff_smooth = single(nan(size(dff)));
for i=1:size(dff,2)
    Fbuffer(:, 1:end-1) = Fbuffer(:, 2:end);
    Fbuffer(:,end) = dff(:,i);
    dff_smooth(:, i) = single(nanmean(Fbuffer, 2));
end

%%
%Inputs: 
%frames_per_reward_range
%cov_bool
reward_per_frame_range  = 1./frames_per_reward_range;
E1_mean                 = mean(analyze_mean(E1_sel));
E1_std                  = sqrt((E1_sel/num_E1)'*analyze_cov*(E1_sel/num_E1));
E2_subord_mean          = zeros(num_E2,1);
E2_subord_std           = zeros(num_E2,1); 
E1_analyze              = dff(:,E1_sel); 
E2_analyze              = dff(:,E2_sel); 
for E2_i = 1:num_E2
    subord_sel                      = E2_sel;
    subord_sel(E2_sel_idxs(E2_i))   = 0; 
    E2_subord_mean(E2_i)            = mean(analyze_mean(subord_sel));     
    var_i                           = subord_sel'*analyze_cov*subord_sel; 
    E2_subord_std(E2_i)             = sqrt(var_i);     
end

E2_sum_analyze = sum(E2_analyze,2); 

%signals needed for target detection:
cursor_obs                      = dff_smooth*decoder; 
E1_meadff                 = mean(E1_analyze,2);
E2_meadff                 = mean(E2_analyze, 2); 
E1_mean_max                     = max(E1_meadff); 
[E2_dom_samples, E2_dom_sel]    = max(E2_analyze, [], 2);
E2_subord_meadff          = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);


%%
%Iterate on T value, until perc correct value is achieved using truncated
%neural activity

%T:
min_prctile         = tset.cb.E2mE1_prctile; %A good default is 98
T0                  = max(cursor_obs);
T                   = T0; 
T_min               = prctile(cursor_obs, min_prctile);

%E2:
E2_coeff0           = 0.5;
E2_coeff            = E2_coeff0; %multiplies the std dev, for figuring out E2_subord_thresh
E2_coeff_min        = 0.05; 
E2_subord_thresh    = E2_subord_mean+E2_subord_std*E2_coeff;

%E1:
E1_coeff0           = 0;
E1_coeff            = E1_coeff0;
E1_thresh           = E1_mean + E1_coeff*E1_std; %E1_mean_max; %E1_mean;


T_delta = 0.05;
E2_coeff_delta = 0.05; %0.05 
E1_coeff_delta = 0.05; %0.05 
task_complete = 0;

T_vec = []; 
E2_coeff_vec = []; 
E1_coeff_vec = []; 

reward_per_frame_vec = []; 

max_iter = 10000;
iter = 0;

%%
%If using data covariance:
rand_num_samples = 500000;
while(~task_complete)
    T_vec           = [T_vec T];
    E2_coeff_vec    = [E2_coeff_vec E2_coeff]; 
    E1_coeff_vec    = [E1_coeff_vec E1_coeff];    

    
    %1) E2-E1 > alpha
    c1 = find(cursor_obs >= T); 
    %2) E1 < mu
    c2 = find(E1_meadff <= E1_thresh);
    %3) E2_subord > mu (anded with previous constraint)
    %For each idx, subtract the 
    c3 = find(E2_subord_meadff >= E2_subord_thresh(E2_dom_sel)); 
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
    reward_prob_per_frame   = sum(hits_valid)/length(dff);    
    %----------------------------------------------------------------------
    
    reward_per_frame_vec = [reward_per_frame_vec reward_prob_per_frame]; 
   
    %Update T:
    if((reward_prob_per_frame >= reward_per_frame_range(1)) && (reward_prob_per_frame <= reward_per_frame_range(2)))
        task_complete = 1;
        disp('target calibration complete!');
    elseif(reward_prob_per_frame > reward_per_frame_range(2))
        %Task too easy, make T harder:
        T = T+T_delta; 
    elseif(reward_prob_per_frame < reward_per_frame_range(1))
        %Task too hard, make T easier:
        T = T-T_delta; 
        %If we swept the full range of T, lower E2_coeff, reset T:
%         if(T<T_min)
%             disp('T through full range, lower E2_coeff, reset T:'); 
%             T=T0; 
%             E2_coeff = E2_coeff - E2_coeff_delta;
%             E2_subord_thresh = E2_subord_mean+E2_subord_std*E2_coeff; 
%         end
%         %If we swept the full range of E2, increase E1_coeff, reset E2:
%         if(E2_coeff < E2_coeff_min)
%             disp('E2 coeff through full range, increase E1_coeff, reset E2_coeff:'); 
%             E2_coeff = E2_coeff0;
%             E1_coeff = E1_coeff + E1_coeff_delta;
%             E1_thresh = E1_mean + E1_coeff*E1_std; %E1_mean_max; %E1_mean;
%         end
    end
%     T
%     E2_coeff
%     E2_subord_thresh
    iter = iter+1;
    if(iter == max_iter)
        task_complete = 1;
        disp('Max Iter reached, check reward rate / baseline data'); 
    end
end

%%
%Summary results of cal: 
disp(['T: ', num2str(T)]); 

num_hits_no_b2base = length(hit_idxs_no_b2base);
disp(['num baseline hits WITHOUT B2BASE: ' num2str(num_hits_no_b2base)]); 

num_valid_hits = length(valid_hit_idxs);
disp(['num valid hits (WITH B2BASE): ' num2str(num_valid_hits)]); 


% b2base_thresh = 0.5*T;
% hits_valid = ones(length(hit_times),1); 
% if length(hit_times) > 1
%     for i = 2:length(hit_times)
%         b2base_bool = sum(cursor_obs(hit_times(i-1):hit_times(i)) <= b2base_thresh) > 1; 
%         hits_valid(i) = b2base_bool; 
%     end
% end
% disp('num baseline hits WITH B2BASE:'); 
% num_valid_hits = sum(hits_valid)

% valid_hit_times = hit_times(find(hits_valid)); 


cursor_amp = (max(cursor_obs)-min(cursor_obs));
cursor_offset = cursor_amp/10; 
max_cursor = max(cursor_obs); 

%%
% Calculate parameters for auditory feedback
[fb_cal]        = define_fb_calibration(cursor_obs, fbset, T);

%
%%
%Save the results: 
%1) All the steps here
%2) Just the target parameters for running BMI

%1)All the steps here
clear h
save_path = fullfile(savePath, ['target_calibration_' var_str '.mat']); 
target_cal_ALL_path = save_path; 
save(save_path, 'num_hits_no_b2base', 'num_valid_hits'); 

%2)Just the target parameters for running BMI
target_info_file = ['BMI_target_info_' var_str '.mat'];
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

function [fb_cal] = define_fb_calibration(cursor_obs, fbset, T)

%Assumes cursor = E2-E1, and cursor_target is positive.
%Maps cursor to auditory feedback.
% freq = a*exp(b*(cursor_trunc-cursor_min))
fb_cal.settings = fbset;

%Calculate:
fb_cal.cursor_min         = prctile(cursor_obs, fbset.min_perctile); 
%min(cursor_obs); %for fb, cursor is ceil to this value
fb_cal.cursor_max         = T; %for fb, cursor is floor to this value
fb_cal.cursor_range       = fb_cal.cursor_max - fb_cal.cursor_min; 
% % freq = a*exp(b*(cursor_trunc-cursor_min))
fb_cal.a                  = fb_cal.settings.freq_min; 
fb_cal.b                  = (log(fb_cal.settings.freq_max) - log(fb_cal.a))/fb_cal.cursor_range; 
end
