function [tset] = define_BMI_task_settings()

%Imaging environment file for baseline acquisition
tset.baseline_env = ...
    fullfile('G:\VivekNuria\utils', 'Tseries_baseline_15.env');
%Imaging environment file for BMI acquisition
tset.bmi_env = ...
    fullfile('G:\VivekNuria\utils', 'Tseries_BMI_40.env');


%
%% imaging
tset.im.frameRate                    = 29.752; 
tset.im.zoom                         = 1.5; 
tset.im.posz                         = 0;
tset.im.chan_data                    = struct('label', 'g', ...
                                                          'chan_idx', 2); %in RGB, G is 2nd

%% rois:
% Parameters for auto cell detection:
% Following were for zoom=2 on bruker soma:
% template_diam = 25; %diamter of difference of Gaussians in pixels
% thres = 0.5; %cell detection threshold as correlation coefficient
% cell_diam = 7; %CELL_DIAM is diameter used for dilation.
% finemode = 1; %imTemplateMatch will be used instead of normxcorr2. It will be slower.
% temmode = 0;  % 0 is for full circle (soma) 1 is for donuts (membrane)
% Following were for zoom=1.5 on bruker soma
tset.roi.template_diam              = 15; %diamter of difference of Gaussians in pixels
tset.roi.thres                      = 0.6; %cell detection threshold as correlation coefficient
tset.roi.cell_diam                  = 14; %CELL_DIAM is diameter used for dilation.
tset.roi.finemode                   = 1; %imTemplateMatch will be used instead of normxcorr2. It will be slower.
tset.roi.temmode                    = 1;  % 0 is for full circle (soma) 1 is for donuts (membrane)

%% calibration:  
tset.cb.sec_per_reward_range   = [100 70]; 
tset.cb.baseline_len           = 7.5*60; %seconds
tset.cb.f0_win_bool            = 1; %during cb, 
%estimate f0 using the window
tset.cb.dff_win_bool           = 1;
tset.cb.f0_init_slide          = 0; %during cb, 
%if 0, f0 is only used after f0_win samples.  if 1, f0 is
%adapted in the window from 0 to f0_win samples.
tset.cb.E2mE1_prctile          = 98;  


%% bmi: 
tset.bmi_len                   = 40*60; %seconds
tset.prefix_win                = 100; 
%set this to 'nonBufferUpdateCounter', 'initFrameBase', number samples to ignore at start of bmi acqn

tset.f0_win                    = 2*60*ceil(tset.im.frameRate); 
%Period at the beginning without BMI to establish BL 

tset.dff_win                   = 10; 
%'movingAverageFrames', number of frames to use for smoothing dff

tset.range_norm_bool           = 1; 
%normalize each neuron by its range

tset.cursor_zscore_bool        = 0; 
%- if 1, neural activity is zscored before going into
%cursor calculation. if 0, neural activity is not zscored.  

tset.rewardDelayFrames         = 10; 
%TODO confirm arduino code triggers reward immediately, so it doesn't add
%extra delay

tset.back2BaseFrameThresh      = 2; %2 frames of back2base 
tset.relaxationTime            = 0; 
tset.b2base_coeff              = 0.5; 
%a frame counts as back2base if cursor < b2base_coeff*T, where T is target cursor value.  

%% random stim

tset.rs.IHSImean = 15; 
tset.rs.IHSIrange = 15; 

%%
%Baseline inputs: 
%baseline_frameRate

%%
%Baseline outputs: 
%frames_per_reward_range
%

