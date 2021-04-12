function [RX_files] = my_get_and_save_rx_hm(RX_Array,...
    RX_AziAngle,RX_EleAngle, phiFinal, thetaFinal, TX_RX_AziAngle, ...
    TX_real_pos_AziAngle, pat_azel_conv,el_coarse, az_shift, MEAS, tempDataFolder)
%UNTITLED Summary of this function goes here
%   pat_azel_conv linear pattern (non-dB)

path = tempDataFolder;
if ~exist(path, 'dir')
    mkdir(path);
end

cmin = -92;
cmax = -65;
ymin = -30;
ymax = 60;
th = -92;

RX_files = {};

for ii=1:length(RX_Array)
    [RSSFinal] = get_sim_hm_for_meas(RX_AziAngle,RX_EleAngle, RX_Array(ii), pat_azel_conv,el_coarse, phiFinal, az_shift);
    %RSSFinal = [RSSFinal, RSSFinal(:,1)];
    % make_angle_heatmap(phiFinal, el_coarse, RSSFinal, cmin, cmax, ymin, ymax, th)
    % title(sprintf('Sim RX %s for TX azimuth %f', MEAS, TX_real_pos_AziAngle(ii) - TX_real_pos_AziAngle(1)))
    filename = ['HM_' MEAS '_TXPOS' num2str(ii) '.mat'];
    finalFileName = fullfile(path, filename);
    RX_files{ii} = finalFileName;
    save(finalFileName, 'phiFinal', 'thetaFinal', 'RSSFinal');
end



end

