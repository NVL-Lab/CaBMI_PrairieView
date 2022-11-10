function [im_s, min_val_return, max_val_return] = scale_im(im, min_perc, max_perc)

min_val = prctile(im(:), min_perc); 
im_s = im-min_val; 
im_s(im_s < 0) = 0; 
max_val = prctile(im_s(:), max_perc); 
im_s(im_s>max_val) = max_val; 
im_s = double(im_s)/double(max_val); 

min_val_return = min_val;
max_val_return = max_val + min_val; 