function [img] = my_imagesc_xscaled(az,el,data,clims,xmin, xmax, th)
%imagesc can only work with equidistant azimuth angles. This one can work
%with any azimuth angles
%   Detailed explanation goes here

az_comp = xmin:1:xmax;

data_comp = zeros(size(el,2),size(az_comp,2));

% fill up the new data with 1 degree granularity
% choose the closest one for each field
for azidx=1:size(az_comp,2)
    azdiff = abs(az-az_comp(azidx));
    [cmin, closest] = find(azdiff == min(azdiff));
    if size(closest,2)>1 || size(closest, 1)>1
        warning('Imagesc_xscaled found two closest neighbors');
        closest = closest(1);
        cmin = cmin(1);
    end
    data_comp(:,azidx) = data(:,closest);
end


img = imagesc(az_comp, el, data_comp, clims);
img.AlphaData = double(logical(~isnan(data_comp) - (data_comp<th)));

end

