function [] = write_mat2csv_frank(TxFolderLists,RxPosLists, ...
    generateSingleFile, generateAggregateFile,...
    writeSinglePathPrefix, writeAggregatePathPrefix,...
    tempDataFolder ,dataType, varargin)
    % write .mat to .csv
    %
    %Input param:
    %   TxFolderLists:  list of directories in which the Tx .dat is saved,  will include 'O'
    %   RxPosLists: n * 2 matrix, indicating the position of the rx
    %   generateSingleFile: boolean, used to specify whether or not will generate single file
    %   generateAggregateFile: boolean, used to specify whether or not will generate aggregate file,
    %   which means that for each rx position, will aggregate all possible tx power available for that rx
    %   writeSinglePathPrefix: path prefix for the single files
    %   writeAggregatePathPrefix: path prefix for aggregate file
    %   tempDataFolder: path fot the temporary data folder
    %   writeSinglePathPrefix: path prefix to which the data will be saved
    %   dataType: str, either 'raytracing' or 'meas'
    %   
    
%     tempName = string.empty;
%     for i = 1 : length(TxFolderLists)
%         tempName(i) = convertCharsToStrings(TxFolderLists(i));
%     TxFolderLists = tempName;

    if (generateSingleFile == false && generateAggregateFile == false)
        error("at least one boolean generate* must be true!");
    end

    if (generateSingleFile == true)
        % check if writeSinglePathPrefix already exists, if not create one
        if ~exist(writeSinglePathPrefix, 'dir')
            mkdir(writeSinglePathPrefix);
        end
    end

    if (generateAggregateFile == true)
        if ~exist(writeAggregatePathPrefix, 'dir')
            mkdir(writeAggregatePathPrefix);
        end
    end

    
    % load antenna pattern based on whether raytracing or measurement
    if strcmp(dataType, 'raytracing')
        load('SiversIMATRX_Matrix.mat');

        % load('FM25240_15.mat');
        % % need to add additional 
        % pat_azel = pat;        
        % [az, el, resmat] = getData(meas_tuple);
    elseif strcmp(dataType, 'meas')
        load('FM25240_15.mat');
        % need to add additional 
        pat_azel = pat;
        % % pick the first experiment data
        % meta = dataRaw.TX_RX_pair.meas_runs(1);
        % selected_meas = meta.get_beam_type_measurements('pencil', 'pencil');
        % [az, el, resmat] = getData(selected_meas);
        % filename = strcat(filename, '_meas')
    end
    
    [m,n] = size(RxPosLists);
    for counter = 1 : m
        firstVisited = true;
        for filePathIndex = 1 : length(TxFolderLists)
            TX_pos = TxFolderLists(filePathIndex);
            load(fullfile(TX_pos,'0p1_new.mat'));
            load(fullfile(TX_pos,'TxRx_Angles.mat'));

            [parentDir,tx_pos,ext] = fileparts(TX_pos);
            


            RX_pos_x = RxPosLists(counter,1);
            RX_pos_y = RxPosLists(counter,2);
            
            RX_pos = Receiver_Ray(RX_pos_x,RX_pos_y);
            
            RX_AziAngle = Rx_AziAngle(RX_pos_x,RX_pos_y,:);
            RX_AziAngle=nonzeros(RX_AziAngle)';
            [rx_r, rx_c] = size(RX_AziAngle);
            if (rx_c == 0)
                fprintf("TX_pos %s, RX_pos (%d, %d) generates null RX_AziAngle, pass\n", tx_pos, RX_pos_x, RX_pos_y);
                continue
            end
            RX_EleAngle = Rx_EleAngle(RX_pos_x,RX_pos_y,:);
            RX_EleAngle=nonzeros(RX_EleAngle)';
            
            TX_RX_AziAngle = Tx_AziAngle(RX_pos_x,RX_pos_y,:);
            TX_RX_AziAngle=nonzeros(TX_RX_AziAngle)';
            
            TX_RX_EleAngle = Tx_EleAngle(RX_pos_x,RX_pos_y,:);
            TX_RX_EleAngle=nonzeros(TX_RX_EleAngle)';
            
            % load azimuth-elevelation antenna angle ranges 
            load('AzEl_Angles.mat');




            %% Definitions for plotting

            cmin = -92;
            cmax = -65;
            ymin = -60;
            ymax = 60;
            th = -100;


            %% Apply the radiation pattern on the TX side

            % The raytracer found the paths to be at TX_RX_AziAngle and TX_RX_EleAngle
            % In this script, tx angles are often used in the range 0 to 180
            tx_angles = tx_angles + 90;

            TX_real_pos_aziangles = TX_RX_AziAngle;
            % We plot the RX heatmaps assuming a fixed TX elevation angle of -15°
            % TX_real_pos_eleangles = -15 * ones(length(TX_RX_A_EleAngle),1);

            % or round to best link opportunity
            TX_real_pos_eleangles = round(TX_RX_EleAngle(1)) * ones(length(TX_RX_EleAngle),1);

            % Apply the transmitter radiation pattern
            [RX_Array] = apply_tx_rad_pat(TX_RX_AziAngle,TX_RX_EleAngle, ...
                TX_real_pos_aziangles, TX_real_pos_eleangles, RX_pos, pat_azel);


            % %% RX - Get heatmaps for fixed TX angles
            rx_pos = sprintf('RX_%d_%d', RX_pos_x, RX_pos_y);

            % scale to roughly -180 to 180
            phiFinal = az_coarse - min(az_coarse(:))-180+1;
            thetaFinal = el_coarse;
            if (generateAggregateFile == true)
                if (filePathIndex == 1 || firstVisited == true)
                    phiOld = phiFinal;
                    thetaOld = thetaFinal;

                elseif (~isequal(phiOld, phiFinal) || ~isequal(thetaOld, thetaFinal))
                    error("the phiOld or thetaOld does not match with the current!\n");
                end
            end
            % horizontal shift to center the heatmap to LOS
            % use this to find the offset to the LOS direction
            az_shift = 0;

            % Get the RX heatmaps for each fixed TX angle
            [RX_files] = my_get_and_save_rx_hm(RX_Array,...
                RX_AziAngle,RX_EleAngle, phiFinal, thetaFinal, TX_RX_AziAngle, ...
                TX_real_pos_aziangles, 10.^(pat_azel/10),thetaFinal, az_shift, rx_pos, tempDataFolder);

            %% RX - Get Max RSS heatmap over all TX angles
            % get max RSS across all heatmaps

            [RSSmax, TXRSSmaxAngles, theta, phi] = ...
                my_max_RSS_per_RX_angle(RX_files, TX_RX_AziAngle - TX_RX_AziAngle(1), el_coarse, th);

            phiFinal = phi;
            thetaFinal = theta;
            RSSFinal = RSSmax;

            if (generateAggregateFile == true)
                if (filePathIndex == 1 || firstVisited == true)
                    % aggregate RSS
                    RSS_agg = RSSFinal;
                    firstVisited = false;
                else
                    [row, col] = size(RSS_agg);
                    for r = 1 : row
                        for c = 1 : col
                            if (~isnan(RSS_agg(r,c)))
                                if (~isnan(RSSFinal(r,c)))
                                    % add the power
                                    % convert to mW, add then turn to dBm
                                    RSS_agg(r,c) = 10 * log10(10^(RSS_agg(r,c)/10) + 10^(RSSFinal(r,c)/10));
                                end
                            elseif(~isnan(RSSFinal(r,c)))
                                % choose the value from RSSFinal
                                RSS_agg(r,c) = RSSFinal(r,c);
                            end
                        end
                    end
                end
            end
            % TX_RX_angles = TX_RX_AziAngle - TX_RX_AziAngle(1);

            % write single files
            if (generateSingleFile == true)
                filename = strcat('TX_', tx_pos ,'_' ,rx_pos ,'_' ,dataType);
                
                % write files
                thetaSuffix = "thetaFinal.csv"
                phiSuffix = "phiFinal.csv"
                RSSSuffix = "RSSFinal.csv"
                parentDir = fullfile(writeSinglePathPrefix, filename)
        
                thetaFile = fullfile(writeSinglePathPrefix, filename, thetaSuffix);
                phiFile = fullfile(writeSinglePathPrefix, filename, phiSuffix);
                RSSFile = fullfile(writeSinglePathPrefix, filename, RSSSuffix);
                if ~exist(parentDir, 'dir')
                    mkdir(parentDir)
                end
                writematrix(thetaFinal, thetaFile, 'Delimiter', 'tab');
                writematrix(phiFinal, phiFile, 'Delimiter', 'tab' );
                writematrix(RSSFinal, RSSFile, 'Delimiter', 'tab');
                fprintf('write %s.mat to directory : %s\n', filename, parentDir);
            end
        end

        if (generateAggregateFile == true)
            filename = strcat(rx_pos ,'_' ,dataType);
            
            % write files
            thetaSuffix = "thetaFinal.csv"
            phiSuffix = "phiFinal.csv"
            RSSSuffix = "RSSFinal.csv"
            parentDir = fullfile(writeAggregatePathPrefix, filename)
    
            thetaFile = fullfile(writeAggregatePathPrefix, filename, thetaSuffix);
            phiFile = fullfile(writeAggregatePathPrefix, filename, phiSuffix);
            RSSFile = fullfile(writeAggregatePathPrefix, filename, RSSSuffix);
            if ~exist(parentDir, 'dir')
                mkdir(parentDir)
            end
            writematrix(thetaOld, thetaFile, 'Delimiter', 'tab');
            writematrix(phiOld, phiFile, 'Delimiter', 'tab' );
            writematrix(RSS_agg, RSSFile, 'Delimiter', 'tab');
            fprintf('write %s.mat to directory : %s\n', filename, parentDir);
        end        
    end
end
    

% function [az, el, resmat]  = getData(meas_tuple)
%     az = sort(unique([meas_tuple.phi_rx]));
%     el = sort(unique([meas_tuple.theta_rx]));
%     RSS = [meas_tuple.rss];

%     %% Max RSS (Data Acquisition)
%     resmat = NaN*zeros(length(az), length(el));
%     for ii = 1:length(RSS)
%         pidx = az == meas_tuple.phi_rx(ii);
%         tidx = el == meas_tuple.theta_rx(ii);
%         if isnan(resmat(pidx, tidx))
%             % guess multiple times there are same theta and phi, and choose the
%             % largest rss among them
%             % pidx, tidx each will only have exactly one place which is 1,
%             % resmat(pidx, tidx) will select that posbill
%             resmat(pidx, tidx) = RSS(ii);
%         elseif RSS(ii) > resmat(pidx, tidx)
%             resmat(pidx, tidx) = RSS(ii);
%         end
%     end    
% end
