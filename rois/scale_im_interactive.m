function [im_sc_struct, num_im_sc] = scale_im_interactive(im, im_sc_struct, num_im_sc)
%INPUT:
%OUTPUT: 
screen_size = get(0,'ScreenSize');
%%
%%Prompt user to change min_perc, max_perc
scale_complete_bool = 0; 
param_vec = [0 100]; 
while(~scale_complete_bool)
    disp('Current [min_perc max_perc]:')	 
    disp(param_vec); 
    param_vec = input('enter [min_perc max_perc]: ');
    while(length(param_vec) ~= 2)
        disp('Error in input!'); 
        param_vec = input('enter [ min_perc max_perc]: ');
    end
    min_perc = param_vec(1); 
    max_perc = param_vec(2); 
    [im_sc, im_min, im_max] = scale_im(im, min_perc, max_perc);    
    
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(im_sc), colormap bone
    axis square
    title(['min prc: ' num2str(min_perc) '; max prc: ' num2str(max_perc)]); 
    
    %----------------------------------------------------------------------
    in = input('Want to Store this Scaling? y/n:   ', 's');
    in = lower(in);
    if(isempty(in) || strcmp(in, 'y'))    
        %Check you don't already have it: 
        already_added = 0; 
        for i = 1:num_im_sc
            if(sum(im_sc_struct(i).minmax ~= param_vec) ==0)
                already_added = 1; 
                disp('already saved it!'); 
                break;
            end
        end
        
        if(~already_added)
            num_im_sc = num_im_sc+1; 
            im_sc_struct(num_im_sc).im = im_sc; 
            im_sc_struct(num_im_sc).minmax_perc = param_vec; 
            im_sc_struct(num_im_sc).minmax = [im_min im_max]; 
            im_sc_struct(num_im_sc).min = im_min; 
            im_sc_struct(num_im_sc).max = im_max; 
            im_sc_struct(num_im_sc).min_perc = min_perc; 
            im_sc_struct(num_im_sc).max_perc = max_perc;       
        end
    end
    
    in = input('Want to Continue More Scalings? y/n:   ', 's');
    in = lower(in);    
    if(strcmp(in, 'n'))
        scale_complete_bool = 1;
    end    
end