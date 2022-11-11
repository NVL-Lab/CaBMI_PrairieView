function autoscale(im, hAx)

if nargin < 2
    hAx = imgca;
end

if nargin < 1
    im = getimage(hAx);
end

im = im;
im = im(isfinite(im)); %excludes NaN, Inf

if isa(im, 'double') | isa(im, 'single')
    immax = max(im(:));
    immin = min(im(:));
    im = (im - immin) / (immax - immin);
end

im = im(~isnan(im));
tolerance = [0.005];
sl = stretchlim(im, tolerance);
range = getrangefromclass(im);
clim = sl * range(2);
%clim(1) = 0;
if isa(im, 'double') | isa(im, 'single') 
   clim = clim * (immax - immin) + immin;
end
set(hAx, 'CLim', clim);