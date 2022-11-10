function [roi_data] = roi_data2chan(roi_data)
% num_chan = 2; 
num_chan = length(roi_data.chan);
for chan_i = 1:num_chan
    roi_data.chan(chan_i).idxs      = find(roi_data.chan_logical(chan_i,:)); 
    roi_data.chan(chan_i).num_rois  = length(roi_data.chan(chan_i).idxs); 
    %Build roi_mask, roi_mask_bin
    roi_data.chan(chan_i).roi_mask      = zeros(roi_data.num_rows,roi_data.num_cols); 
    roi_data.chan(chan_i).roi_mask_bin  = zeros(roi_data.num_rows,roi_data.num_cols); 
    roi_data.chan(chan_i).im_roi        = roi_data.im_bg; 
    for sel_i = 1:length(roi_data.chan(chan_i).idxs)
        i_sel = roi_data.chan(chan_i).idxs(sel_i); 
        roi_sel     = roi_data.roi_bin_cell{i_sel}; 
        roi_idxs    = find(roi_sel);
        
        roi_data.chan(chan_i).roi_mask(roi_idxs) = i_sel; 
        roi_data.chan(chan_i).roi_mask_bin(roi_idxs) = 1; 
        roi_data.chan(chan_i).im_roi(:,:,3) = roi_data.chan(chan_i).roi_mask_bin;
    end
end