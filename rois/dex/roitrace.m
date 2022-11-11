function rtrace = roitrace(mat, binroi)

sm = size(mat); sm(end+1:3) = 1;
mat = reshape(mat, [sm(1)*sm(2) sm(3:end)]);

binroiflat = reshape(binroi, [size(binroi,1)*size(binroi,2) size(binroi,3)]);

for i = 1:size(binroi, 3)
    roiindices = find(binroiflat(:,i));
    tomean = mat(roiindices, :,:,:,:,:,:,:,:);
    if numel(roiindices) == 0, tomean = NaN; end
    if any(isnan(tomean))
        rtrace(i,1,:,:,:,:,:,:,:) = nanmean(tomean, 1);
    else
        rtrace(i,1,:,:,:,:,:,:,:) = mean(tomean, 1);
    end
end

%keyboard
%size(tomean)

