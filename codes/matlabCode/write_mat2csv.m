function [] = write_mat2csv(matFileLists, writePathPrefix,dataType, varargin)
% write .mat to .csv
%
%Input param:
%   matFileLists:  list of paths
%   writePathPrefix: path prefix to which the data will be saved
%   dataType: str, either 'raytracing' or 'meas'
%   

% check if writePathPrefix already exists, if not create one
if ~exist(writePathPrefix, 'dir')
    mkdir(writePathPrefix);
end

for filePathIndex = 1 : length(matFileLists)
        filePath = matFileLists(filePathIndex);
        dataRaw = load(filePath);
        [filepath,filename,ext] = fileparts(filePath);
        if strcmp(dataType, 'raytracing')
            [az, el, resmat] = getData(meas_tuple);
        elseif strcmp(dataType, 'meas')
            % pick the first experiment data
            meta = dataRaw.TX_RX_pair.meas_runs(1);
            selected_meas = meta.get_beam_type_measurements('pencil', 'pencil');
            [az, el, resmat] = getData(selected_meas);
            filename = strcat(filename, '_meas')
        end
        % write files
        thetaSuffix = "thetaFinal.csv"
        phiSuffix = "phiFinal.csv"
        RSSSuffix = "RSSFinal.csv"
        parentDir = fullfile(writePathPrefix, filename)

        thetaFile = fullfile(writePathPrefix, filename, thetaSuffix);
        phiFile = fullfile(writePathPrefix, filename, phiSuffix);
        RSSFile = fullfile(writePathPrefix, filename, RSSSuffix);
        if ~exist(parentDir, 'dir')
            mkdir(parentDir)
        end
        writematrix(el, thetaFile, 'Delimiter', 'tab');
        writematrix(az, phiFile, 'Delimiter', 'tab' );
        writematrix(resmat', RSSFile, 'Delimiter', 'tab');
        fprintf('write %s.mat to directory : %s\n', filename, parentDir);
    end


function [az, el, resmat]  = getData(meas_tuple)
    az = sort(unique([meas_tuple.phi_rx]));
    el = sort(unique([meas_tuple.theta_rx]));
    RSS = [meas_tuple.rss];

    %% Max RSS (Data Acquisition)
    resmat = NaN*zeros(length(az), length(el));
    for ii = 1:length(RSS)
        pidx = az == meas_tuple.phi_rx(ii);
        tidx = el == meas_tuple.theta_rx(ii);
        if isnan(resmat(pidx, tidx))
            % guess multiple times there are same theta and phi, and choose the
            % largest rss among them
            % pidx, tidx each will only have exactly one place which is 1,
            % resmat(pidx, tidx) will select that posbill
            resmat(pidx, tidx) = RSS(ii);
        elseif RSS(ii) > resmat(pidx, tidx)
            resmat(pidx, tidx) = RSS(ii);
        end
    end    
end

end