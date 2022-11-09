%%%TODO'S!!!!
%{
Nuria:
- define the BNCs hookups --> Frame Trigger
checking something
Base: disp number of target hits (without b2base)
BMI: disp number of delivered holo stims
streamline analysis of baseline, pretrain, bmi data, plots in same format
rewards per minute plot
anticipatory licking?
Baseline CA vs stim time
 Close-loop E2 not being high-ish
analyze baseline E2
 what window around stim triggers reward
%
%For Experiment:
%A protocol for determing reliably stimmed cells, and determing power and
duration of stim
%}

%% Main protocol for the experiment
%--------------------------------------------------------------------------
%BEFORE ANIMAL IN BOX:
%DO:
% Hook up BNCs: 
% 1) BMI solenoid, AI5
% 2) Monaco Trig, AI6
% 3) Frame Trig, AI7
% 4) Holo Trig PFI1
% 
% Power Arduino: 
%   (Power supply needed to power solenoid, can't control solenoid on USB power)
% Voltage Recording: All Inputs Active (check 6+7)
%
%Fill syringe with sucrose cuz of gravity
%Run 'main_test_190923.m'
%   check nidaq pulses are received in voltage_rec: 1) frame trigger, 2) trig photostim (monaco), 3)
%   trig reward (bmi solenoid / arduino)
%   calibrate solenoid opening time
%
%load pyctrl expt for the mouse
%In load cell: remove offset. Collect sensor baseline data
% convert to csv
% put the file in pre folder
%
%Put mouse in
%put gel from headbar to ear
%adjust spout so mouse can lick
%put objective
%--------------------------------------------------------------------------

%% DEFINE PATHS
%--------------------------------------------------------------------------
%DO:
%Input 'folder', as directory to write to.
%--------------------------------------------------------------------------
[task_settings] = define_BMI_task_settings();
[fb_settings]   = define_fb_audio_settings();
fb_bool = 0; %If any experiment will need audio
%Initialize arduino:
if(fb_bool) 
    a = arduino(fb_settings.arduino.com, ...
        fb_settings.arduino.label);
else
    a = [];
end

%     a = arduino(task_settings.fb.arduino.com, ...
%         task_settings.fb.arduino.label);

cd G:\VivekNuria\Code\HoloBMI
%DEFINE PATH_DATA: 
%
%LOAD PATHS: 
load_path = define_and_load_bmi_paths();

%SAVE PATHS: 
home_dir = 'G:\VivekNuria\Code\HoloBMI';
cd(home_dir)
env_dir = 'G:\VivekNuria\utils';

% define Animal, day and folder where to save
animal = 'NVI17'; day = 'D19';
folder = 'E:\holobmi_E\191123';
savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
path_data.load_path = load_path; 
path_data.home_dir = home_dir; %home_dir
path_data.env_dir = env_dir; %contains env files for prairie
path_data.savePath = savePath; 
path_data.im = fullfile(savePath, 'im'); %directory for imaging data
if ~exist(path_data.im, 'dir')
    mkdir(path_data.im);
end

connectivity_bool = 0;

%%
%DO: enter zoom (either 1.5 or 2)
zoom = 2; 
posz = 0;
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
frameRate = 29.989;
chan_data = struct(...
    'label', 'g', ...
    'chan_idx', 2); %in RGB, G is 2nd
num_chan = length(chan_data); 
onacid_bool = false;


%%
%{
%--------------------------------------------------------------------------
%DO:
%SUMMARY: 
Option 1:
Turn up imaging power, turn on average of 32 frames, save the image
Option2:
%Take 1000 frame video, calc mean, stddev images, use for
%selecting ROI's.
%
%find FOV
%-disable the motor control!!!!
%Option 1 instructiosn: TODO
%Option 2 instructions:
%take 1000 frame video
%convert it using image-block ripping utility
%load the converted tif into imagej
% [mean, std dev image]: Image->Stack->ZProject, choose mean, choose std
% dev.  
% Put the path into:
%   im_summary_path

%TODO: make an acquisition for this.
%--------------------------------------------------------------------------
%}
%% Select + Scale summary image (to choose ROI)
%--------------------------------------------------------------------------
%DO:
%-Input 'im_summary_path'
%-Run 'scale_im_interactive' : choose percentiles for scaling image
%--------------------------------------------------------------------------

option1_bool = 1; 
if option1_bool
    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    disp('Connecting to prairie')
    pause(2);    
    im_summary = pl.GetImage_2(2, px, py);
    pl.Disconnect();
else
    im_summary_path    = ...
        fullfile('G:\vivek\190822_NY35_good_stim_tests\NY35\D1_test', 'green_std.tif'); 
    % fullfile('E:\vivek\190822\NY35\D1_test', 'chan_mean.tif'); 
    exist(im_summary_path)
    im_summary = imread(im_summary_path); 
end

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

%% REMEMBER TO TURN OFF PHASE OFFSET
%% TURN OFF THE MANIPULATOR 
%% TURN OFF AUTOSCALE

%%
%--------------------------------------------------------------------------
%DO:
%-Input index to 'im_sc_struct'
%--------------------------------------------------------------------------
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

%% INIT ROI_DATA
%--------------------------------------------------------------------------
%DO:
%-Confirm parameters for automatically identifying ROI
%-Run this cell
%-if template matching sucks, can set 'auto_init' to 0
%--------------------------------------------------------------------------
auto_init = 1;  %initializes roi_data using automatic cell detection: 
% Parameters for auto cell detection:
% Following were for zoom=2 on bruker soma:
% template_diam = 25; %diamter of difference of Gaussians in pixels
% thres = 0.5; %cell detection threshold as correlation coefficient
% cell_diam = 7; %CELL_DIAM is diameter used for dilation.
% finemode = 1; %imTemplateMatch will be used instead of normxcorr2. It will be slower.
% temmode = 0;  % 0 is for full circle (soma) 1 is for donuts (membrane)
template_diam = 15; %diamter of difference of Gaussians in pixels
thres = 0.6; %cell detection threshold as correlation coefficient
cell_diam = 14; %CELL_DIAM is diameter used for dilation.
finemode = 1; %imTemplateMatch will be used instead of normxcorr2. It will be slower.
temmode = 1;  % 0 is for full circle (soma) 1 is for donuts (membrane)
if auto_init
    %FIND ROI AUTOMATICALLY 
    [mask_intermediate, ~] = imFindCellsTM (im_bg, template_diam, thres, cell_diam, finemode, temmode);
    init_roi_mask = bwlabel(mask_intermediate);
    findCenter (init_roi_mask, im_bg);
    roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, chan_data);
else
    roi_data = init_roi_data(im_bg, num_chan, chan_data);
end

% use the time now to draw some neurons in the screen
% if result looks weird check pixel size, zoom and parameters of the
% function
% if rois are too small, remove them and draw them yourself

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
% roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, chan_data);
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
    %HoloMask
    holoMask = roi_data.roi_mask; 
    screen_size = get(0,'ScreenSize');
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);    
    imagesc(holoMask)
    axis square; 
    title(['holoMask num roi: ' num2str(roi_data.num_rois)]); 
    
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
roi_mask = roi_data.roi_mask;
save(roi_data_file, 'roi_mask', 'plot_images', 'im_sc_struct', 'roi_data'); 

% filetosave = fullfile(savePath, 'red.mat');
% load(filetosave)
% h = figure;
% imshow(holoMask)

%% prepare SEQUENTIAL HOLO STIM of individual neurons
%Load environment to prairie
% load environment
env_path = fullfile(env_dir, "Tseries_VivekNuria_holo_all.env")
exist(env_path)
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + env_path;
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

%% %Create gpl, xml files for individual points
%DEFINE GPL FILE
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
init_markpoints = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', 100, ...
    'SpiralSize', 0.7, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals

markpoints_data = repmat(init_markpoints, [roi_data.num_rois 1]); 
for roi_i = 1:roi_data.num_rois
    spiral_size_um  = 2*roi_data.r(roi_i)*micronsPerPixel.x;
    spiral_size     = spiral_size_conversion*spiral_size_um;
    markpoints_data(roi_i).SpiralSize = spiral_size; %double(ceil(10*spiral_size)/10);
end
% creates holos
createGplFile_v2(savePath, markpoints_data, roi_data.x, roi_data.y, posz, roi_data.r, px, zoom)

%% XML: Sequential Single Cell Stim
xml_seq_path = fullfile(savePath, 'seq_single_stim.xml'); 
%{
Params Summary:
-num_sequences

-Initial Delay: time (ms) between each stimulation
-power
-numSpirals
-Repetitions
-Iter - number of iterations
-IterDelay - time beteen iterations
-InterPointDelay
-

%}

numberNeurons= roi_data.num_rois; %max(max(holoMask));
% numberNeurons=max(max(holoMask));
num_sequences = 1;
stim_sequence = repmat(1:numberNeurons, [1 num_sequences]); 

power_conversion = 0.004; %0.2 -> 50, 0.4->100

seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
%
InitialDelay = 2000; %(ms) time bw stim delivery
seq_stim_params.InitialDelayVector = InitialDelay*ones(1,numberNeurons);
%
power = 40;
power_converted = power*power_conversion;
seq_stim_params.PowerVector = power_converted*ones(1,numberNeurons);
%2
Duration = 30;
seq_stim_params.DurationVector = Duration*ones(1,numberNeurons);
% 
numSpirals = 10;
seq_stim_params.SpiralVector = numSpirals*ones(1,numberNeurons);
%
Repetitions = 1;
seq_stim_params.RepetitionsVector = Repetitions*ones(1,numberNeurons);
%Darcy sometimes recommends increasing 'Repetitions' and decreasing
%Duration.  This changes the distribution of the spirals in time over the
%cell.
%
InterPointDelay =  0.12;
seq_stim_params.InterPointDelayVector = InterPointDelay*ones(1,numberNeurons); 

% seq_stim_params
createXmlFile_sequential_single_cell(xml_seq_path, seq_stim_params, stim_sequence);

% Update prairie view repetitions based on num neurons to stim
stim_time_per_neuron = InitialDelay/1000+InterPointDelay;
num_reps_seq_stim = ceil(numberNeurons*stim_time_per_neuron*frameRate);
len_seq_stim = numberNeurons*stim_time_per_neuron/60;

disp(['Number of Repetitions in PrairieView: ' num2str(num_reps_seq_stim)])
disp(['Stim time per neuron (s): ' num2str(stim_time_per_neuron)]); 
disp(['Num neurons: ' num2str(numberNeurons)]); 
disp(['Length (min): ' num2str(len_seq_stim)])


% if position of stim cells looks different, "smaller/bigger" check the
% pixel size

%--------------------------------------------------------------------------
%DO: 
% upload .gpl in MarkPoints (Top half)
% upload .xml in MarkPoints (Bot half)
% update T-series repetitions in Prairie View with above number
% Make sure Voltage Recording has all channels enabled
% Make sure you turn on the laser power and pmt's
% while running paint the neurons
% TODO: automate uploading
%--------------------------------------------------------------------------

%% Run HOLO STIM to check stim-able neurons
%This stims one neuron at a time.
%--------------------------------------------------------------------------
%DO: 
%1) Do live scan check
%--------------------------------------------------------------------------

clear s
expt_str = 'holostim_seq'; %previously 'holostim' 
mask = roi_data.roi_mask;
expectedLengthExperiment = ceil(num_reps_seq_stim*1.5); 
HoloAcqnvsPrairie_v2(path_data, expt_str, mask, expectedLengthExperiment)
% HoloAcqnvsPrairie(folder, animal, day, holoMask)
%TODO: make this closed loop, and wait for the neurons to be inactive
%before stimming them.
%{
%--------------------------------------------------------------------------
%DO: 
%Image-block ripping utility: Convert the holostim acqn (2 files)
%Load holostim_seqX.mat
%Load voltage recording for plotting.  
%   Use 'Import Data' in matlab
%   Import csv to matlab: output type is Numeric Matrix.  
%   Name it "voltageRec")
%--------------------------------------------------------------------------
%}
%%
%(Image-Block Ripping Utility) Convert holostim file with bruker converter 
% load the VoltageRec to check the results of holoStim
min_duration = 40; %stims can't occur within this number of samples of voltageRec
plot_win = 1000; 
plotHoloStimTimeLock(holoActivity, voltageRec, min_duration, plot_win)
%ToDo: allow us to select the idxs of neurons to plot
%ToDo: for plotting, do sliding window deltaf/
%ToDo: a version that just plots each individual neuron, we type 'Y' or 'N'
%to make it a candidate


%%
%--------------------------------------------------------------------------
%DO: 
%Summary: choose good candidate stim neurons based on single stim.
%Enter E2_base
% (I often choose more than 4 neurons, manually stim the neurons.
% then re-run once you've chosen your 4.)
%--------------------------------------------------------------------------
E2_candidate =[11 6 31 28]; %  46 3 29 8 13 17 15 
% E2_base = sort([21    36   127   196], 'ascend')

%% Holo stim of Ensemble neurons
% Make GPL (points), BOT (measure activity)
% -select markpoints_data fo
close all
sel_idxs = unique(E2_candidate); 
[sel_roi_data, sel_idxs] = select_roi_data(roi_data, sel_idxs); 
sel_markpoints_data = markpoints_data(sel_idxs); 

%--
%GPL for Stim Ensemble
% savePrairieFiles(savePath, pl, 'GPL_candidates_')
% creates holos
gpl_candidates_path = createGplFile_v2(savePath, sel_markpoints_data, sel_roi_data.x, sel_roi_data.y, posz, sel_roi_data.r, px, zoom, 'GPL_candidates_');
% createGplFile(savePath, StimMask, posz, px, 'ensemble_')
%BOT
bot_candidates_path = fullfile(savePath, 'BOT_candidates.cfg'); 
createBot_v2(bot_candidates_path, sel_roi_data.x, sel_roi_data.y, sel_roi_data.r)
% createBot(savePath, x(E2_base),y(E2_base))
%NOTE: 
%If it can't be loaded in prairie because prairie says the file is in use
%by another program, it can be because matlab didn't release the file.
%Re run, saving to a different filenaem

%XML:
numberNeurons = length(E2_candidate);
ens_default_params.numberNeurons = numberNeurons;
%SPECIFY:
ens_default_params.PowerVector          = 30*power_conversion*ones(1,numberNeurons);
ens_default_params.DurationVector       = 20*ones(1,numberNeurons);
ens_default_params.RepetitionsVector    = 1*ones(1,numberNeurons);
ens_default_params.SpiralVector         = 10*ones(1,numberNeurons);
%Don't Change:
ens_default_params.InitialDelayVector   = 0*ones(1,numberNeurons);
ens_default_params.UncagingLaser        = "Monaco"; 
ens_default_params.AllPointsAtOnce      = "False"
ens_default_params.Iter                 = 1; %how many times to go through and stim each cell.
ens_default_params.IterDelay            = 0; %Time (ms) between iterations
InterPointDelay                         =  0.12;
ens_default_params.InterPointDelayVector = InterPointDelay*ones(1,numberNeurons);
xml_candidates_path = fullfile(savePath, 'XML_candidates.xml'); 
createXmlFile_basic(xml_candidates_path, ens_default_params);

pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile(path_data.env_dir, "Tseries_VivekNuria_holo_4.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect();

% Check neuron size, sometimes they are genourmous sometimes they are tiny.
% Improvise, adapt, overcome

%--------------------------------------------------------------------------
%DO: 

%1) upload the GPL file
%2) upload XML file
%3) Click BOT in Image Window, load BOT.cfg, only display ROI of interest
%4) Run BOT, and adjust (duration, power, repetitions) of each neuron's stim
%5) Note down below the params for each neuron
%--------------------------------------------------------------------------

%%
%{
%DO:
%In the GUI Mark Points Series, remove bad neurons 
%   1) from Points/Groups, 
    2) from Mark Point Series
    3) from BOT 
%For chosen neurons, Enter Stim Parameters you like in GUI.
EXPORT each of the above, with the following names: 
%}

%%
%--------------------------------------------------------------------------
%DO: 
%1) Choose E2_base, if different than E2_candidate
%2) Specify: PowerVector, DurationVector, RepetitionsVector, SpiralVector 
% else:

gpl_path = fullfile(savePath, 'GPL_ens.gpl');
xml_path = fullfile(savePath, 'XML_ens.xml');
bot_path = fullfile(savePath, 'BOT_ens.cfg');

copyfile(gpl_candidates_path, gpl_path)
copyfile(bot_candidates_path, bot_path)
copyfile(xml_candidates_path, xml_path)


%--------------------------------------------------------------------------
% % E2_base = E2_base([2 11 9 15])
E2_base = E2_candidate; 
% %GPL (define the points):
% sel_idxs = unique(E2_base); 
% [stim_roi_data, stim_idxs] = select_roi_data(roi_data, sel_idxs);
% stim_markpoints_data = markpoints_data(stim_idxs); 
% StimMask = stim_roi_data.roi_mask;
% numberNeurons = length(E2_base);
% createGplFile_v2(savePath, stim_markpoints_data, ...
%     stim_roi_data.x, stim_roi_data.y, posz, stim_roi_data.r, px, 'GPL_ens_')
% 
% %XML:
% ens_stim_params.numberNeurons = length(E2_base);
% %SPECIFY:
% ens_stim_params.PowerVector = [15 35]*power_conversion;
% ens_stim_params.DurationVector = [20 20];
% ens_stim_params.RepetitionsVector = 1*[1 1];
% ens_stim_params.SpiralVector = [10 10];
% %Don't Change:
% ens_stim_params.InitialDelayVector = 0*ones(1,numberNeurons);
% ens_stim_params.UncagingLaser = "Monaco"; 
% ens_stim_params.AllPointsAtOnce = "False"
% ens_stim_params.Iter = 1; %how many times to go through and stim each cell.
% ens_stim_params.IterDelay = 0; %Time (ms) between iterations
% InterPointDelay =  0.12;
% ens_stim_params.InterPointDelayVector = InterPointDelay*ones(1,numberNeurons);
% xml_ens_path = fullfile(savePath, 'XML_ens.xml'); 
% createXmlFile_basic(xml_ens_path, ens_stim_params);
% 
% %BOT
% botPath = fullfile(savePath, 'BOT_ens.cfg'); 
% createBot_v2(botPath, stim_roi_data.x, stim_roi_data.y, stim_roi_data.r)




%% Baseline acquisition
%Note: loads the result of OnAcid / holoMask
%Do this after we confirm we can stim some cells
%--------------------------------------------------------------------------
%DO: 
%Remove Red Channel from Image Window 1 (prairie view).
%0) (zero pmt+power) put water
% check FOV didn't move
close all
imshow(im_bg)
%1) start video
%2) start load cells
%3) run and start pyctrl
%4) Run following cell

%0) if you need to convert other files from previous experiments do it now
% before starting baseline check that the SNR is good an keep a screenshot
% of it

%% --------------------------------------------------------------------------
if ~onacid_bool
    AComp = 0;
else
    load(roi_data_file);
end
% Baseline environment already removes MARKPOINTS and set the reps to 27000
holoMask = roi_data.roi_mask;
[base_mat_path, base_dat_path] = ...
    BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMask, task_settings);
% BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMaskRedGreen, onacid_bool, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activityim_bg

%--------------------------------------------------------------------------
%D0:
%0) Abort T-series (cuz of voltage recording)
%1) pyctrl stop
%2) load cells stop
%3) video stop
%4) Drag load cell data to folder
%5) Drag video to folder
%--------------------------------------------------------------------------

%% Selection of neurons
% plots neurons so we can select which ones we like the most 

%Copy paste base_file path: 
base_file = base_mat_path; 
%Can manually enter a previous path: 
% base_file = fullfile(savePath, 'BaselineOnline190526T113422.mat')
if onacid_bool
    totalneurons = min(size(AComp,2), 20);
else
    totalneurons = max(max(roi_data.num_rois));
    CComp = [];
    YrA = []; 
end
load(base_file); 
% totalneurons = 40; 
% plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)
plotNeuronsBaseline(baseActivity, CComp, YrA, 30)
%TODO:  
%ToDo: for plotting, do sliding window deltaf/f
%%
%--------------------------------------------------------------------------
%D0:
%1) Choose E1_base
%--------------------------------------------------------------------------
%when choosing E1s, avoid neurons that are bright or next to bright ones

% 
%Manually enter and confirm the BMI neurons:
% E2_candidate = unique([9 15 23 29]); %unique also sorts
% E2_base = sort([7 9 12 29]; %8 21 10 6 17 3 18
E1_base = sort([25 15 10 26], 'ascend')  % 34 17 7 28
ensembleNeurons = [E1_base, E2_base];
plotNeuronsEnsemble(baseActivity, ensembleNeurons, [ones(1,length(E1_base)) 2*ones(1,length(E2_base))])
select_roi_data(roi_data, [E2_base, unique(E1_base)]); 
% E2_candidates = [39 45 59 37 88 6 26 46 78 48 22 20 33]
%E2_candidate =[19 33 37 24]; % 25 16 32 15  36 35 18 12  
%% for E3 experiments:
%E2_candidate =[8 19 24 12]; %  19 13 4 16 2 34 15 12 9 14

% E3_base = unique([8 19 24 12]); % 22 21 16 23 4 5 9 
% E2_base = sort([22 13 4 16], 'ascend'); % 
% % 


%%
%OPTION: Use previously collected BMI data as the baseline data: 
%
% bmi_file = fullfile(savePath, 'BMI_online190515T010526.mat'); 
% bmi_data = load(bmi_file); 
% bmi_base = fullfile(savePath, ['base_' 'BMI_online190515T010526.mat']);
% baseActivity = bmi_data.data.bmiAct(:, ~isnan(bmi_data.data.bmiAct(1,:))); 
% save(bmi_base, 'baseActivity'); 
%
% E1_base = [1 2 3 4]; 
% E2_base = [5 6 7 8]; 

%% Calibrate Target with Baseline simulation
%--------------------------------------------------------------------------
%D0: (nothing)
%1) Parameters: 
% - sec_per_reward_range
% - f0_win (F0: how many frames to average over)
% - dff_win (F for Dff: how many frames to average over)
%--------------------------------------------------------------------------

% base_file = fullfile(savePath, 'BaselineOnline190514T221822.mat')
% base_file = bmi_base; 

exist(base_file)
n_f_file = base_file;
ndata = load(n_f_file);
num_base_samples = sum(~isnan(ndata.baseActivity(1,:))); 
baseline_frameRate = num_base_samples/(15*60);
A_file = roi_data_file; %fullfile(savePath, 'red.mat'); 
exist(A_file)
onacid_bool = 0

sec_per_reward_range = [120 90]; 
% sec_per_reward_range = [10 5]


frames_per_reward_range = sec_per_reward_range*baseline_frameRate;
disp('Time (s) per reward range: '); 
disp(sec_per_reward_range); 
disp('Frames per reward range: '); 
disp(frames_per_reward_range)
% sec_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 

E2mE1_prctile = 98; 
target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 2*60*ceil(frameRate)
dff_win_bool = 1
dff_win = 4
 
reward_per_frame_range = 1./frames_per_reward_range

cursor_zscore_bool = 0;
f0_init_slide = 0; 

close all
[target_info_path, target_cal_ALL_path, fb_cal] = baseline2target_vE1strict_fb(n_f_file, A_file, onacid_bool,  ...
    E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
    prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
    cursor_zscore_bool, f0_init_slide, E2mE1_prctile, fb_settings);

% [target_info_path, target_cal_ALL_path] = baseline2target_vE1strict(n_f_file, A_file, onacid_bool,  ...
%     E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
%     prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
%     cursor_zscore_bool, f0_init_slide, E2mE1_prctile);

%ToDo: return the filename
% run the simulation of baseline
%To Do: Show the percent correct of the pretrain period, based on the
%calibration. 


% CAREFUL!!! You don't want T too small (all noise) or too big (stim can not reach it)
% ideal value is   0.2 < T < 0.5
% if T too low, I can only recommend to change the E1 neurons or E2
% if T too high you can change the ensemble neurons or this E2mE1_prctile = 95;
% sometimes the E2mE1 oercentil is no good if the E2 neurons are very
% active. that is why is better to choose, active but not crazy active
% neurons
%--------------------------------------------------------------------------
%D0:
%Note down: 
% - T value
% T = 0.31
% num_valid_hits:7
% num_hits: 91
%--------------------------------------------------------------------------
%% Holo stim checking connectivity
% create randomize run for each individual neuron of the ensemple
%--------------------------------------------------------------------------
%D0:
%-Manually enter: powerVector, durationVector
%--------------------------------------------------------------------------

%% Runs connectivity
%--------------------------------------------------------------------------
%D0:
% choose 'num_conn' = number of times to stim each neuron
%--------------------------------------------------------------------------
if connectivity_bool
    if pl.Connected()
        pl.Disconnect();
    end
    savePrairieFiles(savePath, pl, 'connectivity_pre')
    numberNeurons = stim_roi_data.num_rois;
    num_conn = 10;
    conn_sequence = repmat(1:numberNeurons, [1 num_conn]);  
    %Randomize order:
    conn_sequence = conn_sequence(randperm(length(conn_sequence))); 
    
    ens_conn_params = ens_stim_params; 
    ens_conn_params.IterDelay = 3000; %Time (ms) between iterations    
    InitialDelay = 3000; %(ms) time bw stim delivery
    ens_conn_params.InitialDelayVector = InitialDelay*ones(1,numberNeurons);    
    
    % seq_stim_params
    xml_conn_pre_path = fullfile(savePath, 'conn_pre.xml');     
    createXmlFile_sequential_single_cell(xml_conn_pre_path, ens_conn_params, conn_sequence);

    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    loadCommand = "-tsl " + fullfile(path_data.env_dir, "Tseries_VivekNuria_holo_4.env");
    pl.SendScriptCommands(loadCommand);
    pl.Disconnect();
end

%% Run 
%--------------------------------------------------------------------------
%D0:
%1) (zero PMT, power) Add water to imaging window if needed
%2) Load xml
%--------------------------------------------------------------------------
if connectivity_bool
    clear s
    StimMask = stim_roi_data.roi_mask; %reload just in case
    expectedLengthExperiment = 3*7000;
    ConnectivityAcqnvsPrairie_v2(savePath, expectedLengthExperiment, StimMask, 'PRE')
end

%% create stims for pretrain
%DO:
%Change iterations to 121 in the mark point series

% iterations              = 121;
% pretrain_xml_path       = fullfile(path_data.savePath, 'XML_pretrain.xml'); 
% pretrain_params         = ens_stim_params; 
% pretrain_params.Iter    = iterations; 
% createXmlFile_basic(pretrain_xml_path, pretrain_params);
% createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector, iterations, 'preTrain', false)

%%
%Compute vectorHolo
%--------------------------------------------------------------------------
%D0:
%1) Confirm IHSI mean, range
%2) seedBase - if we will seed the baseline, then set to 1.  
% - if seedBase 0, we wait for baseline before starting stims
%--------------------------------------------------------------------------

frameRate = 30 %baseline_frameRate
baseFrames = 2*60*30; 
expectedLengthExperiment = 70*60*frameRate

% IHSImean, IHSIrange
IHSImean = 20; 
IHSIrange = 10; 
[vectorHolo, ISI] = createVectorHolo(frameRate, expectedLengthExperiment, IHSImean, IHSIrange, false);

seedBase = 0; %Set this to 1 if you will seed the baseline
if ~seedBase
    vectorHolo = vectorHolo + baseFrames;
end
% num imaging reps should be 75600 = 72000+3600

%% Load Ensemble BOT:  create masks bot and image to check during experiment
%In BOT, load 'BOT_ens.cfg'
%If it was some how modified:
make_ens_bot = 0; 
if make_ens_bot 
    BOT_ens_path = fullfile(path_data.savePath, 'BOT_ens.cfg'); 
    createBot_v2(BOT_ens_path, stim_roi_data.x,stim_roi_data.y,stim_roi_data.r);
end

%% run Pre-training

%%
%Seed BMI baseVal using Pretrain
%--------------------------------------------------------------------------
%D0: (only seeded)
%1) seedBase - if we will seed the baseline, then set to 1. 
% - if seedBase 0, we wait for baseline before starting stims
%2) Copy-paste BMI_target_info filename (into 'pretrain_file')
%--------------------------------------------------------------------------
seedBase = 0; 
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
if seedBase
    %TODO:
    pretrain_file = 'BMI_online190523T010653'
    load(fullfile(savePath, pretrain_file)); 
    pretrain_base = data.baseVector; 
    pretrain_base(:, isnan(pretrain_base(1,:))) = [];
    baseValSeed = pretrain_base(:,end)
end

%%
% Pre-training

%--------------------------------------------------------------------------
%D0:
%
%Confirm 'target_info_path'
% %Change the Mark Points:
% %Clear Point Series, Load pretrain.xml
%Make 121 iterations, 

% IMPORTANT put "Wait for Trigger" = First Reptition, Trigger
%Trigger Selection: Start with External, PFI1
%

% Then, before running cell:
%0) zero pmt, laser.  put water under objective.
%1) start video
%2) start load cells
%3) start pyctrl
%--------------------------------------------------------------------------
close all
imshow(im_bg)
clear s
baselineCalibrationFile = target_info_path;
vectorVTA = []
%expt_str: 
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

expt_str = 'HoloVTA_pretrain'; 
debug_bool = 0; 
debug_input = []; 

BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input, baseValSeed);

% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
%     expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input, baseValSeed)

% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v3(folder, animal, day, ...
%     expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input);
%TODO: add seed functionality
%
% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v2(path_data.savePath, path_data.env_dir, ...
%     expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input, baseValSeed);
%saves filename with expt_str
% BMIAcqnvsPrairienoTrialsHoloCL(folder, animal, day, expt_str, baselineCalibrationFile, baseline_frameRate, vectorHolo, vectorVTA, cursor_zscore_bool)

%--------------------------------------------------------------------------
%D0:
%Stop:
%1) pyctrl
%2) load cells
%3) video
%--------------------------------------------------------------------------

%% run BMI
%--------------------------------------------------------------------------
%D0:
% Remove mark points!!!!
% Then, before running cell:
%0) put water under objective
%1) start video
%2) start load cells
%3) start pyctrl
%--------------------------------------------------------------------------


% %Get baseValSeed from HoloVTA_pretrain!  load file, take the last valid
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
task_settings.fb.arduino.duration = 1
playTone(a,...
    task_settings.fb.arduino.pin,...
    fb_freq_i,...
    task_settings.fb.arduino.duration)


%%
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
baselineCalibrationFile = target_info_path;
%%
close all
imshow(im_bg)
%  baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
vectorHolo = [];
vectorVTA= []; 
debug_bool = 0; 
debug_input = []; 
cursor_zscore_bool = 0; 

%T1 = 0.27

expt_str = 'BMI'; 
BMIAcqnvsPrairienoTrialsHoloCL_fb_debug_enable_v4(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input, baseValSeed, fb_bool, fb_cal, a);

% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
%     expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input, baseValSeed)
% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v3(folder, animal, day, ...
%     'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input);


% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v2(path_data.savePath, path_data.env_dir, ...
%     'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input, baseValSeed)

%Stop:
%1) pyctrl
%2) load cells
%3) video

%%
%DO: random reward
%Clear vectorVTA
%Clear vectorHolo
bmi_no_reward_bool = 1; 

if bmi_no_reward_bool

    close all
    imshow(im_bg)
     baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
    vectorHolo = [];
    vectorVTA= []; 
    debug_bool = 0; 
    debug_input = []; 
    cursor_zscore_bool = 0; 

    %T1 = 0.27

    expt_str = 'BMI_no_reward'; 
    BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
        expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
        cursor_zscore_bool, debug_bool, debug_input, baseValSeed)
    % BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v3(folder, animal, day, ...
    %     'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    %     cursor_zscore_bool, debug_bool, debug_input);


    % BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v2(path_data.savePath, path_data.env_dir, ...
    %     'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    %     cursor_zscore_bool, debug_bool, debug_input, baseValSeed)

    %Stop:
    %1) pyctrl
    %2) load cells
    %3) video
end

%%
%DO: random reward
%Define vectorVTA
%Clear vectorHolo
random_reward_bool = 1; 

if random_reward_bool

    close all
    imshow(im_bg)
    %  baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
%     vectorHolo = [];
%     vectorVTA= []; 
    debug_bool = 0; 
    debug_input = []; 
    cursor_zscore_bool = 0; 

    %T1 = 0.27

    expt_str = 'VTA_pretrain'; 
    BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
        expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
        cursor_zscore_bool, debug_bool, debug_input, baseValSeed)
    % BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v3(folder, animal, day, ...
    %     'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    %     cursor_zscore_bool, debug_bool, debug_input);


    % BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v2(path_data.savePath, path_data.env_dir, ...
    %     'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    %     cursor_zscore_bool, debug_bool, debug_input, baseValSeed)

    %Stop:
    %1) pyctrl
    %2) load cells
    %3) video
end


%% HOLO PRETRAIN

expectedLengthExperiment = 40*60*frameRate
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan

IHSImean, IHSIrange
IHSImean = 20; 
IHSIrange = 10; 
[vectorHolo, ISI] = createVectorHolo(frameRate, expectedLengthExperiment, IHSImean, IHSIrange, false);


close all
imshow(im_bg)
clear s
baselineCalibrationFile = target_info_path;
vectorVTA = []
%expt_str: 
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

expt_str = 'Holo_pretrain'; 
debug_bool = 0; 
debug_input = []; 

BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input, baseValSeed)

% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v3(folder, animal, day, ...
%     expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input);
%TODO: add seed functionality
%
% BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v2(path_data.savePath, path_data.env_dir, ...
%     expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input, baseValSeed);
%saves filename with expt_str
% BMIAcqnvsPrairienoTrialsHoloCL(folder, animal, day, expt_str, baselineCalibrationFile, baseline_frameRate, vectorHolo, vectorVTA, cursor_zscore_bool)

%--------------------------------------------------------------------------
%D0:
%Stop:
%1) pyctrl
%2) load cells
%3) video



%% Holo stim checking connectivity
if connectivity_bool
    if pl.Connected()
        pl.Disconnect();
    end
    savePrairieFiles(savePath, pl, 'connectivity_post')
    numberNeurons = stim_roi_data.num_rois;
    num_conn = 10;
    conn_sequence = repmat(1:numberNeurons, [1 num_conn]);  
    %Randomize order:
    conn_sequence = conn_sequence(randperm(length(conn_sequence))); 
    
    ens_conn_params = ens_stim_params; 
    ens_conn_params.IterDelay = 5000; %Time (ms) between iterations    
    InitialDelay = 5000; %(ms) time bw stim delivery
    ens_conn_params.InitialDelayVector = InitialDelay*ones(1,numberNeurons);    
    
    % seq_stim_params
    xml_conn_post_path = fullfile(savePath, 'conn_post.xml');     
    createXmlFile_sequential_single_cell(xml_conn_post_path, ens_conn_params, conn_sequence);

    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    loadCommand = "-tsl " + fullfile(path_data.env_dir, "Tseries_VivekNuria_holo_4.env");
    pl.SendScriptCommands(loadCommand);
    pl.Disconnect();
end

%%
if connectivity_bool
    clear s
    StimMask = stim_roi_data.roi_mask; %reload just in case
    expectedLengthExperiment = 3*7000;
    ConnectivityAcqnvsPrairie_v2(savePath, expectedLengthExperiment, StimMask, 'POST')
end
%% end


%%
%--------------------------------------------------------------------------
%D0:

%1) SAVE THE WORKSPACE IN FOLDER
%2) SAVE THIS Protocol script IN FOLDER (savePath)
%3) Start converting imaging data
%3) Remove mouse
%4) collect load cell sensor baseline data

%--------------------------------------------------------------------------
%%
%NOTES
% rr. During Baseline the laser acted up. Because it was already the 2nd
% baseline and I didn't want to have the animal for another 15min there. I
% did the calibration with the full 15min and taking only the first 10min,
% the results were the same, so I kept the calibration and did NOT redo
% baseline
% second problem, the water went down. I've been refilling it, but it is
% not ideal