function strcMask=obtain_Strc_Mask_from_Mask(mask)
%{
Function to create a structure with a reduce image given a spatial filter 
px,py => shapes of the image
strcMask -> structure with the matrix for spatial filters with px*py*unit
and the positions of that mask
units => index of the neurons in the neuronMask
%}
    
    [x,y]           = find_center(mask);
    roi_ctr         = [x';y']; %size: 2 x num_roi    
    
    %adjust idxs to loop over using 'unique'
    roi_ind = unique(mask(:));
    roi_ind(roi_ind==0) = []; %remove the 0 ind
    num_roi = length(roi_ind);  
    
    strcMask.roi_ind = roi_ind; 
    strcMask.num_roi = num_roi; 
    for u = 1:num_roi
        auxMask = mask;
        auxMask(auxMask ~= roi_ind(u)) = 0;
        posx = find(sum(auxMask,1)~=0);
        posy = find(sum(auxMask,2)~=0);
        strcMask.maxx(u) = posx(end);
        strcMask.minx(u) = posx(1);
        strcMask.maxy(u) = posy(end);
        strcMask.miny(u) = posy(1);
        strcMask.neuronMask{u} = auxMask(posy(1):posy(end), posx(1):posx(end));
        
        %Aux Information:
        strcMask.xctr(u)    = roi_ctr(1,u); 
        strcMask.yctr(u)    = roi_ctr(2,u);
        strcMask.width(u)   = abs(strcMask.maxx(u)-strcMask.minx(u)); 
        strcMask.height(u)  = abs(strcMask.maxy(u)-strcMask.miny(u)); 
    end
end