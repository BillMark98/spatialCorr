function my_make_heat_map(folderNames, folderParentPath, savedFolderPre, suffixes)
    % make heat map based 
%   - folderNames: list of strs, folder where the csv files are saved
%   - folderParentPath: str, the parent dir
%   - savedFolderPre: str, the parent dir of the saved photo
%   - suffixes: list of suffixes for the saved photo file type

%% Definitions for plotting


% default behavior
cmin = -92;
cmax = -65;
ymin = 0;
ymax = 181;
xmin = -180;
xmax = 180;
th = -92;
    % get the last root directory name
    names = split(folderParentPath,"/");
    dirName = names(end);
    % if dir name is specified without an '/' at the end,
    % then the last element is "" choose the next to last element
    if (dirName == "" && length(names) >= 2)
        dirName = names(end - 1);
    end
    thetaFile = "thetaFinal.csv";
    phiFile = "phiFinal.csv";
    RSSFile = "RSSFinal.csv";
    for counter = 1 : length(folderNames)
        folderName = folderNames(counter);
        dataFolder = fullfile(folderParentPath, folderName);

        thetaFolder = fullfile(dataFolder, thetaFile);
        phiFolder = fullfile(dataFolder, phiFile);
        RSSFolder = fullfile(dataFolder, RSSFile);
        if ~isfile(thetaFolder) || ~isfile(phiFolder) || ~isfile(RSSFolder)
            fprintf("the file %s or the two other .csv file does not exist! skip that!\n", thetaFolder);
            continue;
        end

        savedPath = fullfile(savedFolderPre, folderName);
        savedFigureName = strcat(savedPath, ".pdf");
        fprintf("note that in this case already plotted figures wont be plotted!!!!\n");
        if isfile(savedFigureName)
            fprintf("%s already exists, will not be plotted again!\n", savedFigureName);
            continue;
        end

        load(thetaFolder);
        load(phiFolder);
        load(RSSFolder);


        az = phiFinal;
        el = thetaFinal;

        %% Plotting
        fig = figure;
        set(fig, 'Position', [0 0 1000 250])
        th = -100; %Langenfeld=-91, SuperC=-78

        %if sensitivity > th
        %    th = sensitivity;
        %end

        range = max(max(RSSFinal(:))-th, 40);
        img = my_imagesc_xscaled(az, el, RSSFinal, [th th+range], xmin, xmax, th); %-68 for max(RSSFinal,[],'all')
        %img.AlphaData = double(logical(~isnan(RSS) - (RSS<threshold)));
        %imagesc(phiFinal, thetaFinal, rxFinal, clims);
        axis 'image'
        c = colorbar;
        colormap jet
        map=colormap;
        %map(1,:)=1;
        set(gcf, 'Colormap', map)
        set(gca,'YDir','normal')
        yticks([-30 0 30 60])
        ylabel('\theta_{RX} (\circ)')
        ylabel(c, 'RSS (dBm)')
        xlabel('\phi_{RX} (\circ)')
        titleName = sprintf("%s %s", folderName, dirName);
        % replace '_' as ' ' because it will treat the first symbol after _ as a subscript
        % note that to use single quotes!!! not double quotes
        % and do not forget to reassign it to titleName ........
        titleName = strrep(titleName, '_', ' ');
        title(titleName);

        %title('RSS (dBm) depending on angle')

        savedPath = fullfile(savedFolderPre, folderName);
        %% Save

        if ~exist(savedFolderPre, 'dir')
            mkdir(savedFolderPre);
        end
            saveas(gcf, savedPath, 'epsc');
            fprintf("%s saved!\n", savedPath);
            for counter = 1 : length(suffixes)
                saveas(gcf,savedPath,suffixes(counter));
            end
            
            % name = savedPath+".tikz";
            % tikzfile = convertStringsToChars(name);
            % matlab2tikz(tikzfile,  ...
            %     'width', '180pt', 'height', '45pt', ...
            %     'extraaxisoptions', ['unit vector ratio={1 1},'...
            %     'ticklabel style={font={\color{black}\footnotesize}}']);
            % format_tikz_file(tikzfile);
        
            % close all figures
            close all
    end
end