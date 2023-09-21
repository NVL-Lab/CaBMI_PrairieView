
B = roi_data.roi_mask;
G = roi_data.im_bg(:,:,1)/max(max(roi_data.im_bg(:,:,1)));
R = zeros(512,512);
for i=1:4
R(strcMask.yctr(i)-2:strcMask.yctr(i)+2, ...
    strcMask.xctr(i)-2:strcMask.xctr(i)+2)=1;
end
RGB = cat(3,R*255,G*1, B*255);
imshow(RGB)

figure()
B = roi_data.roi_mask/2;
R = zeros(512,512);
G = zeros(512,512);
for i=1:2
R(strcMask.yctr(i)-2:strcMask.yctr(i)+2, ...
    strcMask.xctr(i)-2:strcMask.xctr(i)+2)=1;
end
for i=3:4
G(strcMask.yctr(i)-2:strcMask.yctr(i)+2, ...
    strcMask.xctr(i)-2:strcMask.xctr(i)+2)=1;
end
RGB = cat(3,R,G, B);
imshow(RGB)