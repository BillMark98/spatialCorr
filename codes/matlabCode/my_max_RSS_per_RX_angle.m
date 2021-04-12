function [RSSmax, TXRSSmaxAngles, theta, phi] = ...
    my_max_RSS_per_RX_angle(list_of_files, txangles, RXElAngles, threshold)
%Returns an array with the overall maximum at every position and the
%corresponding TX angles
%   This function loads all the files in the cell array list_of_files one
%   after another. Then, an array that has the dimensions
%   length(RXElAngles)x(360/6+1) is created where each elements represents
%   the measured value for a certain angle combination (elevation,
%   azimuth). The maximum measured value per element is chosen across all
%   files and returned in RSSmax. TXRSSmaxAngles contains the corresponding
%   TX angle for which the maximum was measured. theta and phi are the
%   elevation and azimuth angles as vectors, respectively.


for fileidx=1:length(list_of_files)
    
    load(list_of_files{fileidx})
    
    [RSSFinal, phiFinal, thetaFinal] = ...
    resize_heatmap(RSSFinal, phiFinal, thetaFinal, phiFinal, RXElAngles);
    
    % be compatible with some old files...
    if ~exist('RSSFinal')
        warning('RSS Final does not exist. Use rxUncalibratedFinal instead');
        RSSFinal = rxUncalibratedFinal;
    end

    
    theta = thetaFinal;
    phi = phiFinal;
    
    if exist('RSSFinal')
        if ~exist('RSSmax')
            RSSmax = zeros(size(RSSFinal));
            %RSSmax = zeros(7,61);
            RSSmax(:,:) = NaN;
            TXRSSmaxAngles = zeros(size(RSSFinal));
            TXRSSmaxAngles(:,:) = NaN;
        end
        for tcnt = 1:length(thetaFinal)
            for pcnt = 1:length(phiFinal)
                if (RSSFinal(tcnt, pcnt) > RSSmax(tcnt, pcnt) || ...
                        RSSmax(tcnt, pcnt) == 0 || ...
                        isnan(RSSmax(tcnt, pcnt))) && ...
                        RSSFinal(tcnt, pcnt) > threshold
                    RSSmax(tcnt, pcnt) = RSSFinal(tcnt, pcnt);
                    TXRSSmaxAngles(tcnt, pcnt) = txangles(fileidx);                  
                end
            end
        end
    else
        error('Could not find variable RSSFinal')
    end
    
    clear RSSFinal
        
end


end

