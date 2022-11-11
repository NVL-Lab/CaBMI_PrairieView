%%%TODO'S!!!!
%{
Nuria:
- create the environments for tseries 15 and 40 (for baseline and BMI)
- the CD make it a path from settings
- find the best variables for kenichi in task_settins.image
- check that cleaning saves the baseline
- check sec_per_reward_range is very low (only 1 every 2 min!)
- check that the calibration works well with only C1! and then remove c2
also in dff2cursor_target
- check the Tdelta maybe that's what's up!
- need to clean better baseline2target for now I leave it as is
- check random stim with baseline?
- check data type in computer?
- check first alone, then with a mouse if the stim of d1r fucks up the pmts
- chanching everything to single, check that it doesn't give problems!


in 'define_BMI_tset()'
- define the folder paths!
- check the framerate for what I want
in 'define_and_load_bmi_paths()'
- change the paths
in BMI_Acqnvs_Prairie
- change the paths
- play tone for reward?

%}
%% BEFORE STARTING
%{ 
make sure the BNCs are connected properly
- PCI-AI0 goes to NIDAQ -> P0.1
- PCI-AI1 goes to PFI8 with a T-BNC
the arduino
atm the arduino is sending the activation of the light stim
- digital pin D5 is blue
- digital pin D3 is UV
and the nidaq
the nidaq sends a quick pulse to synchronize each BMI frame
for that we also use voltage recording to get each Prairie-view frame
voltage rec samples/sec min 1000, time(ms) > time_t_series

Need to start the Doric software
- add both channels
- select ttl trigger for both
- amplitude will be: blue (d5):    and UV (d3):
%}

%% DEFINE PATHS and parameters
%--------------------------------------------------------------------------
%DO:
%Input 'folder', as directory to write to.
%--------------------------------------------------------------------------
[tset] = define_BMI_task_settings();
[fbset]   = define_fb_audio_settings();
%Initialize arduino:
if(fbset.fb_bool) 
    a = arduino(fbset.arduino.com, fbset.arduino.label);
else
    a = [];
end

%DEFINE PATH_DATA: 

%SAVE PATHS: 
env_dir = 'E:\Nuria\utils';

% define Animal, day and folder where to save
animal = 'NVI17'; day = 'D19'; date = '221111';
folder = 'E:\D1BMI\';
savePath = fullfile(folder, animal,  date, day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
path_data.env_dir = env_dir; %contains env files for prairie
path_data.savePath = savePath; 
path_data.im = fullfile(savePath, 'im'); %directory for imaging data
if ~exist(path_data.im, 'dir')
    mkdir(path_data.im);
end

%% Get pixel values from prairie

pl = actxserver('PrairieLink.Application');
pl.Connect(); disp('Connecting to prairie...');
% Prairie variables
px = pl.PixelsPerLine();
py = pl.LinesPerFrame();
micronsPerPixel.x = str2double(pl.GetState('micronsPerPixel', 'XAxis')); 
micronsPerPixel.y = str2double(pl.GetState('micronsPerPixel', 'YAxis')); 
pl.Disconnect();
disp('Disconnected from prairie');
% px = 512; 
% py = 512;


num_chan = length(tset.im.chan_data); 


%% Get first image to obtain rois
%{
%--------------------------------------------------------------------------
%DO:
Turn up imaging power, turn on average of 32 frames, save the image
%--------------------------------------------------------------------------
%}

pl = actxserver('PrairieLink.Application');
pl.Connect();
disp('Connecting to prairie')
pause(2);    
im_summary  = pl.GetImage_2(2, px, py);
pl.Disconnect();

%Scale the image, in order to help see ROIs better.
%If no modification to original image needed, just run code, in 
% 'scale_im_interactive' set min_perc = 0, max_perc = 100
im_sc_struct = struct(...
    'im', [], ...
    'minmax_perc', [], ...
    'minmax', [], ...
    'min', [], ...
    'min_perc', [], ...
    'max', [], ...
    'max_perc', []); 
num_im_sc = 0; 
[im_sc_struct, num_im_sc] = scale_im_interactive(im_summary, im_sc_struct, num_im_sc);
close all;

%--------------------------------------------------------------------------
%% INIT ROI_DATA
 
%Can input a different index to choose as the Image for choosing ROI
%Defaults to the last image in 'im_sc_struct'
im_bg = im_sc_struct(end).im; 
h = figure;
imagesc(im_bg); 
axis square
colormap('gray'); 
title('selected background image for identifying ROI'); 
%PLOT_IMAGES data:
%'plot_images' contains a set of images so user can tell if ROI selection is
%appropriate.
plot_images = struct('im', [], 'label', ''); 
plot_images(1).im = im_summary; 
plot_images(1).label = 'green mean';
plot_images(2).im = im_bg; 
plot_images(2).label = 'scaled';

[mask_intermediate, ~] = imFindCellsTM (im_bg, tset.roi.template_diam, tset.roi.thres, ...
    tset.roi.cell_diam, tset.roi.finemode, tset.roi.temmode);
init_roi_mask = bwlabel(mask_intermediate);
find_center (init_roi_mask, im_bg);
roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, tset.im.chan_data);


% use the time now to draw some neurons in the screen
% if result looks weird check pixel size, zoom and parameters of the
% function
% if rois are too small, remove them and draw them yourself
% don't select too many rois!

%%
%Visualize: 
screen_size = get(0,'ScreenSize');
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
hold on;
imagesc(roi_data.im_roi); %colormap('gray');  
axis square;
title('ROI footprint overlay in blue'); 
% scatter(roi_data.x, roi_data.y, pi*roi_data.r.^2, 'r'); 

%
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
% hold on;
imagesc(roi_data.roi_mask); %colormap('gray');  
axis square;
title('ROI Mask'); 
% scatter(roi_data.x, roi_data.y, pi*roi_data.r.^2, 'r'); 

%TESTS
% roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, tset.im.chan_data);
% [roi_ctr] = roi_bin_cell2center_radius(roi_data.roi_bin_cell);

%% Delete ROI if needed
%--------------------------------------------------------------------------
%DO:

%-If auto detected ROI suck, delete ROI
%--------------------------------------------------------------------------
close all;
%Delete ROI if needed: 
disp('Deleting ROIs from image!');
[roi_data] = delete_roi_2chan(plot_images, roi_data);
close all;

%% Add ROI if needed
%--------------------------------------------------------------------------
%DO:
%-Run cell to manually draw additional ROI
%--------------------------------------------------------------------------
disp('Adding ROIs to image!'); 
[roi_data] = draw_roi_g_chan(plot_images, roi_data);
close all

%% SEE ROI if needed
see_roi_data = 1; 
if see_roi_data
    %Roi_Mask
    screen_size = get(0,'ScreenSize');
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);    
    imagesc(roi_data.roi_mask)
    axis square; 
    title(['roi_mask num roi: ' num2str(roi_data.num_rois)]); 
    
    %im_roi:
    screen_size = get(0,'ScreenSize');
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    % hold on;
    imagesc(roi_data.im_roi); %colormap('gray');  
    axis square;
    title(['ROI footprint overlay in blue.  Num ROI: ' num2str(roi_data.num_rois)]);
end

%% Save roi_data

roi_data_file = fullfile(savePath, 'roi_data.mat'); 
save(roi_data_file, 'plot_images', 'im_sc_struct', 'roi_data'); 


%% Baseline acquisition
% Baseline environment already sets the reps to 27000
base_mat_path = Baseline_Acqnvs_Prairie(folder, animal, day, roi_data.roi_mask, tset);
% saves in [savePath, 'baselineActivity'] the activity of all the

%--------------------------------------------------------------------------
%D0:
%0) Abort T-series (cuz of voltage recording)

%% Plot neurons from baseline
% plots neurons so we can select which ones we like the most 

%Copy paste base_file path: 
base_file = base_mat_path; 
%Can manually enter a previous path: 
% base_file = fullfile(savePath, 'BaselineOnline190526T113422.mat')

totalneurons = max(max(roi_data.num_rois));
CComp = [];
YrA = []; 

load(base_file); 
plot_Neurons_Baseline(baseActivity, CComp, YrA, totalneurons)


%% Select MANUALLY ensemble neurons
%Manually enter and confirm the BMI neurons:
E2_base = sort([7 9 12 29], 'ascend'); 
E1_base = sort([25 15 10 26], 'ascend'); 
ensembleNeurons = [E1_base, E2_base];
plot_Neurons_Ensemble(baseActivity, ensembleNeurons, [ones(1,length(E1_base)) 2*ones(1,length(E2_base))])
select_roi_data(roi_data, [unique(E2_base), unique(E1_base)]); 

%%
%OPTION if things went wrong: Use previously collected BMI data as the baseline data: 
%
% bmi_file = fullfile(savePath, 'BMI_online190515T010526.mat'); 
% bmi_data = load(bmi_file); 
% bmi_base = fullfile(savePath, ['base_' 'BMI_online190515T010526.mat']);
% baseActivity = bmi_data.data.bmiAct(:, ~isnan(bmi_data.data.bmiAct(1,:))); 
% save(bmi_base, 'baseActivity'); 
% base_file = bmi_base; 
%
% E1_base = [1 2 3 4]; 
% E2_base = [5 6 7 8]; 

%% Calibrate Target with Baseline simulation
%--------------------------------------------------------------------------
%1) Parameters: 
% - sec_per_reward_range
% - f0_win (F0: how many frames to average over)
% - dff_win (F for Dff: how many frames to average over)
%--------------------------------------------------------------------------
% we can change the basefile again
% base_file = fullfile(savePath, 'BaselineOnline190514T221822.mat')

exist(base_file)
n_f_file = base_file;
ndata = load(n_f_file);
num_base_samples = sum(~isnan(ndata.baseActivity(1,:))); 
baseline_frameRate = num_base_samples/(15*60);
exist(roi_data_file)

sec_per_reward_range = tset.cb.sec_per_reward_range; %[120 90]; 

frames_per_reward_range = sec_per_reward_range*baseline_frameRate;
disp('Time (s) per reward range: '); 
disp(sec_per_reward_range); 
disp('Frames per reward range: '); 
disp(frames_per_reward_range)
% sec_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 


reward_per_frame_range = 1./frames_per_reward_range;

close all
[target_info_path, target_cal_ALL_path, fb_cal] = baseline2target(n_f_file, roi_data_file,  ...
    E1_base, E2_base, frames_per_reward_range, tset, savePath, fbset);


% CAREFUL!!! You don't want T too small (all noise) or too big (stim can not reach it)
% ideal value is   0.2 < T < 0.5
% if T too low, I can only recommend to change the E1 neurons or E2
% if T too high you can change the ensemble neurons or this E2mE1_prctile = 95;
% sometimes the E2mE1 oercentil is no good if the E2 neurons are very
% active. that is why is better to choose, active but not crazy active
% neurons
%--------------------------------------------------------------------------

%%
%Compute vector_stim
%--------------------------------------------------------------------------
%D0:
%1) Confirm IHSI mean, range
%2) seedBase - if we will seed the baseline, then set to 1.  
% - if seedBase 0, we wait for baseline before starting stims
%--------------------------------------------------------------------------

baseFrames = 2*60*30; 
expectedLengthExperiment = 70*60*tset.im.frameRate ;

% IHSImean, IHSIrange 
[vector_stim, ISI] = create_Vector_random_stim(tset.im.frameRate , expectedLengthExperiment, tset.rs.IHSImean, tset.rs.IHSIrange, false);

seedBase = 0; %Set this to 1 if you will seed the baseline
if ~seedBase
    vector_stim = vector_stim + baseFrames;
end
% num imaging reps should be 75600 = 72000+3600

%--------------------------------------------------------------------------

%% run BMI
%--------------------------------------------------------------------------
% DO!!!
% rename the file in the jetball computer!
% OPTIONAL: Get baseValSeed from previous BMI!  load file, take the last valid
% %baseVal
% % load_baseVal = 0; 
% % if load_baseVal
% % baselineCalibrationFile = 'BMI_target_info_20190523T220638.mat';
% pretrain_file = 'BMI_online190524T131817'
% load(fullfile(savePath, pretrain_file)); 
% pretrain_base = data.baseVector; 
% pretrain_base(:, isnan(pretrain_base(1,:))) = [];
% baseValSeed = pretrain_base(:,end)
%% Test FB

fb_freq_i = 7000;
fbset.arduino.duration = 1;
playTone(a,...
    fbset.arduino.pin,...
    fb_freq_i,...
    fbset.arduino.duration)

%%
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan ;
baselineCalibrationFile = target_info_path;
%%
close all
imshow(im_bg)
%  baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan 

% define the type of experiment
expt_str = 'BMI'; % or 'random_stim' or...
BMI_Acqnvs_Prairie(folder, animal, day, ...
    expt_str, baselineCalibrationFile, tset.im.frameRate , vector_stim, [], ...
    tset.cursor_zscore_bool, 0, [], baseValSeed, fb_bool, fb_cal, a);


%--------------------------------------------------------------------------
%D0:

%1) SAVE THE WORKSPACE IN FOLDER
%2) SAVE THIS Protocol script IN FOLDER (savePath)
