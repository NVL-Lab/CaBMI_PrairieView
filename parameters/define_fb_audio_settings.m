function [fb_settings] = define_fb_audio_settings()

fb_settings.fb_bool                = 1; % 1 to have feedback, 0 to silence
fb_settings.target_low_freq        = 1; 
%Set the target cursor value to be the low frequency
fb_settings.freq_min               = 6000; 
fb_settings.freq_max               = 19000; 
fb_settings.arduino.com            = 'COM11';
fb_settings.arduino.label          = 'Mega2560';
fb_settings.arduino.pin            = 'D3';
fb_settings.arduino.duration       = 0.3; %ms, tones update at rate of BMI code, this is the longest a tone will play for
fb_settings.min_perctile            = 90; 


%Added 11.5.19 to have feedback dependent on E1, E2 state
%target_buffer: how much frequency separation there should be between 
%feedback for target achievement vs intermediate feedback
%The following frequencies are arbitrarily decided, could be chosen with principle
if fb_settings.target_low_freq
    fb_settings.target_freq_buffer      = 1000;
    fb_settings.trunc_freq_non_target   = fb_settings.freq_min + fb_settings.target_freq_buffer;
    fb_settings.trunc_freq_E1_state     = 9000; 
else
    fb_settings.target_freq_buffer      = 2000;
    fb_settings.trunc_freq_non_target   = fb_settings.freq_max - fb_settings.target_freq_buffer;
    fb_settings.trunc_freq_E1_state     = 15000; 
end


% added by nuria to solve an error:
fb_settings.lambda_E2mE1 = 0.5; 
fb_settings.lambda_E1 = 0.25; 
fb_settings.lambda_E2 = 0.25;

fb_settings.min_prctile            = 10; %The lowest percentile allowed for E2 minus E1
fb_settings.max_prctile            = 100; %The lowest percentile allowed for E2 minus E1
fb_settings.middle_prctile         = 50; 
fb_settings.obj_max_perctile       = 90; 
