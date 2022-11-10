function freq = cursor2audio(cursor, fb_cal)
%Cursor is the neural control signal.  From old work, this would be: 
%sumE2-sumE1
%
%Auditory feedback calculation
%freq = a*exp(b*(cursor_trunc-cursor_min))
% freq_min = a =  a*exp(0)
% freq_max = a*exp(b*(cursor_max-cursor_min))
% b = log(freq_max/a)/(cursor_max-cursor_min)
%param: cursor_min, cursor_max
%%
%If params.fb.target_low_freq == 1, then we negate cursor and cursor_min,
%cursor_max
%
%This is because cursor is E2-E1.  Our target is gonna be positive.  

%%
%Handle target -> freq:
if fb_cal.settings.target_low_freq == 1
    %This means cursor up makes auditory freq go down:
    cursor      = -cursor; 
    cursor_min  = -fb_cal.cursor_max;
    cursor_max  = -fb_cal.cursor_min;
else
    %This means cursor up makes auditory freq go up:
    cursor_min  = fb_cal.cursor_min;
    cursor_max  = fb_cal.cursor_max;
end

%%
cursor_trunc    = max(cursor, cursor_min); 
cursor_trunc    = min(cursor_trunc, cursor_max); 
freq = fb_cal.a*exp(fb_cal.b*(cursor_trunc-cursor_min));
freq = double(freq); 
% h = figure;
% plot(cursor, freq); 