disp(['Running CrackTeeth2020.m']) ;
%  Does visualization of crack teeth data
%
%  Data are in folder:
%    csv_files
%  from Matt Wilkerson, Oct. 2011
%    1st column is 
%    Remaining Columns are Reads
%    Rows are chromosome location
%
%  Output figures:
%      - Raw data curves
%      - log10 data curves
%      - PCA scatterplot matrix
%      - PCA scatterplot matrix brushed
%      - log10 curves brushed
% Angel Huang Created 9/10/2020

%% Select process to run
ipart = 23;
    %    0 - Read in Raw Data, and Save as .mat file
    %    1 - (In report) Raw Data Curves, original or log scale
    %    2 - (In report) Curve data plot on PC direction
    %    3 - Marginal Distribution plot for each sample (No use)
    %    4 - (In report) PCA scatterplot on 2PC directions, brushed
    %    5 - Projection plot for each sample (No use)
    
    %    2x - Quantile analysis
    %    21 - (In report) Plot raw traces and log scale traces for quantile
    %    22 - (In report) Quantile curve data plot on PC direction
    %    23 - (In report: mean) Marginal Distribution plot of quantile
    %    24 - Quantile PCA scatterplot on 2PC directions, brushed
    
    %    31 - (In report) DWD for crack-size curves
    %    32 - (In report) DiProPerm test on DWD class seperation
    %    33 - (In report) DWD for Quantile representation of data
    %    34 - (In report) DiProPerm test on DWD class seperation of quantiles

%% Read in Raw Data, and Save as .mat file
datSaveName = 'CrackTeeth2020.mat';
%addpath(genpath('C:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\202009_Function data analysis\'));
addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\202009_Function data analysis\'));

if ipart == 0    
    % Navigate to project folder and then execute
    datFolderName = 'synthetic_csvs/';
    fileInfo = dir([datFolderName '*.csv']);
    nFile = numel(fileInfo);    
    dataS = cell(nFile,2);
    dataMat = zeros(948,nFile); % 948 is the max number of clusters
    for iFile = 1:nFile
        splitName = strsplit(fileInfo(iFile).name,'_');        
        teethIDs{iFile,1} = fileInfo(iFile).name(1:end-10);
        teethIDs{iFile,2} = splitName{1};
        teethIDs{iFile,3} = splitName{2};
        teethIDs{iFile,4} = splitName{end-2};
        teethIDs{iFile,5} = splitName{end-1};
        teethIDs{iFile,6} = ~strcmp(teethIDs{iFile,4},'orig');
        crackMask(iFile) = ~strcmp(teethIDs{iFile,4},'orig');
        num = csvread([datFolderName fileInfo(iFile).name]); % data in each file
        dataS{iFile,1} = teethIDs(iFile,1);
        dataS{iFile,2} = num;
        dataS{iFile,3} = log10(num);
        dataS{iFile,4} = num/sum(dataS{iFile,2}); % normalized data (% voxel)
        dataS{iFile,5} = log10(dataS{iFile,4});
        dataMat(1:numel(num),iFile) = num'; % Each column is a sample
    end
    % set up brush matrix
    mcolor = ones(nFile,1) * [0 0 1]; % initialize as blue
    mcolor(crackMask,:) = ones(sum(crackMask),1) * [1 0 0];
    save(datSaveName,'dataS','dataMat','teethIDs','crackMask','mcolor','-v7.3');
end

%%  Load data from previously saved .mat file
load(datSaveName);
cd('./results/');
nFeatures = cellfun(@length,dataS(:,2)); % array of number of features
minFeature = min(nFeatures); % 9
medFeature = median(nFeatures); % 119
maxFeature = max(nFeatures); % 948
legendcellstr = {'cracked','healthy'};
mlegendcolor = [1 0 0; 0 0 1]; % r, b
nTeeth = size(dataS,1);
nCrack = sum(crackMask);
npc = 4;
nHealthy = nTeeth - nCrack;
titlecellstr = {{['n = ' num2str(nTeeth) ' teeth'...
      ' (' num2str(nCrack) ' cracked, ' num2str(nHealthy) ' healthy)']}} ;

% To get some summary
maxFeature_Crack = nFeatures(crackMask);
maxFeature_Healthy = nFeatures(~crackMask);
max(maxFeature_Crack);
min(maxFeature_Crack);
max(maxFeature_Healthy);
min(maxFeature_Healthy);
    
%% Plot raw traces and log scale traces
if ipart == 1
    savestr = [num2str(ipart) '_rawTrace'];
    fig = AH_figure(1,2,savestr);
    colors = 'br'; % corresponding to crackMask = 0,1
    xLim = [0,50];
    for iFile = 1:nTeeth
        subplot(121) % original scale
        l(iFile) = plot(dataS{iFile,2}, colors(crackMask(iFile)+1));
        hold on
        xlim(xLim);ylim([0,800]); xlabel('Feature ID'); ylabel('NumVoxels');
        if iFile == 1; title('Distribution of NumVoxels'); end
        
        subplot(122) % log scale
        plot(log10(dataS{iFile,2}), colors(crackMask(iFile)+1));
        hold on
        xlim(xLim);xlabel('Feature ID'); ylabel('Log_{10}(NumVoxels)');
        if iFile == 1; title('Distribution of Log_{10}(NumVoxels)'); end
    end
    legend([l(1) l(nTeeth)],legendcellstr); 
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
    
    %% Plot histogram of feature number
    savestr = [num2str(ipart) '_nFeatureHist'];    
    H = histogram(nFeatures,20,'FaceColor','r', 'FaceAlpha',0.5);
    % create bin template for both classes
    fig = AH_figure(1,1,savestr);
    histogram(nFeatures(crackMask),'BinEdges',H.BinEdges, 'FaceColor','r', 'FaceAlpha',0.5);
    
    hold on
    histogram(nFeatures(~crackMask),'BinEdges',H.BinEdges, 'FaceColor','b', 'FaceAlpha',0.5);
    
    legend(legendcellstr); 
    title(titlecellstr{:});
    ylabel('numTeeth'); xlabel('numFeatures');
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
    
    %% Plot histogram of feature number
    savestr = [num2str(ipart) '_firstFeatureScatter'];    
    fig = AH_figure(1/2,1,savestr);
    scatter(log(dataMat(1,:)),crackMask);
    ylim([-0.2,1.2]);
    title([titlecellstr{:} 'line = 10000 voxels']);
    vline(log(10000));
    ylabel('isCracked'); xlabel('log(firstFeature)');
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);    
    
end

%% Plot all distributions
if ipart == 2  
    doLog = 1;
    % Prepare data
    nFeature = medFeature;
    mdat = dataMat(1:nFeature,:); % Truncate to shortest sample (min feature count)
    logSuffix = '';
    if doLog == 1
        mdat = log10(mdat+0.1);
        logSuffix = '_log10';
    end    
    % Prepare param
    savestr = [num2str(ipart) '_curvdatPlot_' num2str(nFeature) 'Features' logSuffix] ;
    paramstruct = struct('icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...                         
                         'legendcellstr',{legendcellstr}, ...
                         'mlegendcolor',mlegendcolor, ... % color legend
                         'isubpopkde', 1,...% partition data into subpopulations
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    % Plot all features
    fig = AH_figure(npc,npc,savestr);
    curvdatSM(mdat,paramstruct); % data: dDimension x nSample
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%% Don't need direct marginal distribution plot on variable, but on each case
% % %% Marginal distribution plot
% % % Choose statistics to plot
% % %                      1 - Sample mean
% % %                      2 - Sample Standard Deviation (default)
% % %                      3 - Skewness
% % %                      4 - Kurtosis
% % %                      5 - Median
% % %                      6 - MAD
% % %                      7 - IQR
% % %                      8 - Min
% % %                      9 - Max
% % %                      10 - Range (= max - min)
% % %                      11 - Number of unique values (using Matlab's "unique")
% % %                      12 - Number of most frequent value
% % %                      13 - Number of 0's
% % %                      14 - smallest non-zero spacing
% % %                      15 - Continuity Index (proportion of non-zero 
% % %                                    pairwise distances)
% % %                      16 - Entropy (discrete version)
% % %                      17 - Bowley Skewness (robust version, based on 
% % %                                    quartiles and median)
% % %                      18 - 2nd L-statistic ratio (robust version of variance)
% % %                      19 - 3rd L-statistic ratio (robust version of skewness)
% % %                      20 - 4th L-statistic ratio (robust version of kurtosis)
% % if ipart == 3  
% %     doLog = 1;
% %     % Prepare data
% %     nFeature = minFeature;
% %     mdat = dataMat(1:nFeature,:); % Truncate to shortest sample (min feature count)
% %     logSuffix = '';
% %     if doLog == 1
% %         mdat = log10(mdat+0.1);
% %         logSuffix = '_log10';
% %     end    
% %     % Prepare param    
% %     savestr = [num2str(ipart) '_margDistPlot_' num2str(nFeature) 'Features' logSuffix] ;
% %     nplot = 36;
% %     istat = 9;
% % 
% %     paramstruct = struct('istat',istat,...
% %                          'nplot',nplot,...
% %                          'icolor',mcolor, ...
% %                          'titlecellstr',titlecellstr, ...
% %                          'legendcellstr',{legendcellstr}, ...
% %                          'mlegendcolor',mlegendcolor, ...
% %                          'savestr',savestr, ...
% %                          'iscreenwrite',1) ;
% %     
% %     % Plot marginal distribution plot
% %     fig = AH_figure(npc,npc,savestr);
% %     statstr = MargDistPlotSM(mdat,paramstruct); % data: dDimension x nSample
% %     % Update name of plot with stat performed
% %     savestr = [num2str(ipart) '_margDistPlot_' statstr '_' num2str(nFeature) 'Features' logSuffix] ;
% %     saveas(fig, [savestr '.fig']);
% %     saveas(fig, [savestr '.png']);
% % end

%% Do PCA scatterplot
npc = 4; % number of PCs
doLog = 1;
mcolor = ones(nTeeth,1) * [0 0 1]; % initialize as blue (since scatter plot all data is black too)
mcolor(crackMask,:) = ones(sum(crackMask),1) * [1 0 0];

if ipart == 4
    nFeature = medFeature; % medFeature=119 1000 is the max, can choose smaller value to only consider earlier clusters
    mdat = dataMat(1:nFeature,:); % nFeatures x nSample (each column is a sample)
    mdat(isnan(mdat)) = 0;
    logSuffix = '';
    if doLog == 1
        mdat = log10(mdat+0.1);
        logSuffix = '_log10';
    end
    % Initial PCA
    paramstruct = struct('npc',npc,...
                   'iscreenwrite',1,...
                   'viout',[0 0 0 0 1]);
    outstruct = pcaSM(mdat,paramstruct);
    mpc = getfield(outstruct,'mpc'); % npc x nCluster
                        
    savestr = [num2str(ipart) '_PCAScatPlot_' num2str(nFeature) 'Features' logSuffix] ;    
    paramstruct = struct('npcadiradd',4, ...
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'PC 1'; 'PC 2'; 'PC 3'; 'PC 4'}}, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    fig = AH_figure(npc,npc,savestr);
    scatplotSM(mdat,[],paramstruct);
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%% Marginal distribution plot of each case
if ipart == 5
    doLog = 1;
    doZoom = 1;
    if doLog == 1; logSuffix = '_log10'; end    
    if doZoom == 1
        savestr = [num2str(ipart) '_projPlot_' num2str(nTeeth) 'Samples' logSuffix '_zoom'] ;
    else
        savestr = [num2str(ipart) '_projPlot_' num2str(nTeeth) 'Samples' logSuffix] ;
    end
    %vdir = zeros(zeros(d,1));
    % Plot marginal distribution plot
    fig = AH_figure(4,5,savestr); % number only proximate the size
    %  Make individual distribution plots
    for iTeeth = 1:nTeeth
        % Prepare data
        mdat = dataS{iTeeth,2};
        if doLog == 1
            mdat = log10(mdat+0.1);
        end
        % Prepare param        
        paramstruct = struct('xlabelstr', {['Sample ' num2str(iTeeth)]}) ;
        subplot(5,5,iTeeth);
        projplot1SM(mdat,1,paramstruct); % vdir has the same number of row as mdat
        if crackMask(iTeeth) == 1 % samples are cracked teeth
            set(gca, 'XColor', 'r','YColor','r');
        end
        if doZoom == 1
            ylim([0,8]);xlim([-0.1,1]);
        end
    end
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%% 
if ipart == 6
    doLog = 1;
    nFeature = medFeature; % 1000 is the max, can choose smaller value to only consider earlier clusters
    %vdir = zeros(zeros(d,1));
    % Plot marginal distribution plot
    mdat = dataMat(1:nFeature,:)';
    if doLog == 1
        mdat = log10(mdat+0.1);
        logSuffix = '_log10';
    end 
    savestr = [num2str(ipart) '_featureDistribution_' num2str(nFeature) 'features' logSuffix] ;

    nBin = 40;
    % Get binEdges from all samples
    H1 = histogram(mdat(:,1),nBin, 'FaceColor','r', 'FaceAlpha',0.5);    
    fig = AH_figure(1,3,savestr); % number only proximate the size
    subplot(131)
    histogram(mdat(crackMask,1),'BinEdges',H1.BinEdges, 'FaceColor','r', 'FaceAlpha',0.5);
    hold on
    histogram(mdat(~crackMask,1),'BinEdges',H1.BinEdges, 'FaceColor','k', 'FaceAlpha',0.5);    
    legend(legendcellstr); 
    title('1st feature distribution');
    ylabel('numTeeth'); xlabel('Normalized 1st feature value');
    
    figure();H2 = histogram(mdat(:,2),nBin, 'FaceColor','r', 'FaceAlpha',0.5); 
    set(0,'CurrentFigure',fig)

    subplot(132)
    histogram(mdat(crackMask,2),'BinEdges',H2.BinEdges, 'FaceColor','r', 'FaceAlpha',0.5);
    hold on
    histogram(mdat(~crackMask,2),'BinEdges',H2.BinEdges, 'FaceColor','k', 'FaceAlpha',0.5);    
    title('2nd feature distribution');
    ylabel('numTeeth'); xlabel('Normalized 2nd feature value');
    
    figure();H3 = histogram(mdat(:,1)./mdat(:,2),nBin, 'FaceColor','r', 'FaceAlpha',0.5);
    set(0,'CurrentFigure',fig)
    subplot(133)
    H = histogram(mdat(crackMask,1)./mdat(crackMask,2),'BinEdges',H3.BinEdges, 'FaceColor','r', 'FaceAlpha',0.5);
    hold on
    histogram(mdat(~crackMask,1)./mdat(~crackMask,2),'BinEdges',H3.BinEdges, 'FaceColor','k', 'FaceAlpha',0.5);    
    title('1st/2nd feature distribution');
    ylabel('numTeeth'); xlabel('Normalized 1st/2nd feature value');
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end




%% Calculate quantile
vprob = [0.01:0.01:1]; % 100 intervals
vquant = zeros(numel(vprob),nTeeth);
for iTeeth = 1:nTeeth
    mdat = dataS{iTeeth,2}'; % 1d column vector
    vquant(:,iTeeth) = cquantSM(mdat,vprob,1); % nQuantile x nTeeth
end


%% Plot raw traces and log scale traces
if ipart == 21
    savestr = [num2str(ipart) '_Q_rawTrace'];
    fig = AH_figure(1,2,savestr);
    colors = 'br'; % corresponding to crackMask = 0,1
    xLim = [0,0.95];
    for iFile = 1:nTeeth
        subplot(121) % original scale
        l(iFile) = plot(vprob,vquant(:,iFile)',colors(crackMask(iFile)+1));
        hold on
        
        if iFile == 1; title('Distribution of Quantile(NumVoxels)');
            xlim(xLim); xlabel('Quantile'); ylabel('Quantile(NumVoxels)');
        end
        
        subplot(122) % log scale
        plot(vprob,log10(vquant(:,iFile)'), colors(crackMask(iFile)+1));
        hold on        
        if iFile == 1; title('Distribution of Log_{10}(Quantile(NumVoxels))');
            xlim(xLim); xlabel('Quantile'); ylabel('Log_{10}(Quantile(NumVoxels))');
        end
    end
    legend([l(1) l(nTeeth)],legendcellstr); 
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%% Plot all distributions
if ipart == 22  
    doLog = 1;
    % Prepare data
    nFeature = size(vquant,1); % number of quantiles
    mdat = vquant; % Truncate to shortest sample (min feature count)
    logSuffix = '';
    if doLog == 1
        mdat = log10(mdat+0.1);
        logSuffix = '_log10';
    end    
    % Prepare param
    savestr = [num2str(ipart) '_Q_curvdatPlot_' num2str(nFeature) 'Features' logSuffix] ;
    paramstruct = struct('icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'legendcellstr',{legendcellstr}, ...
                         'mlegendcolor',mlegendcolor, ...
                         'isubpopkde', 1,...% partition data into subpopulations                        
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    % Plot all features
    fig = AH_figure(npc,npc,savestr);
    curvdatSM(mdat,paramstruct); % data: dDimension x nSample
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%% Marginal distribution plot
% Choose statistics to plot
%                      1 - Sample mean
%                      2 - Sample Standard Deviation (default)
%                      3 - Skewness
%                      4 - Kurtosis
%                      5 - Median
%                      6 - MAD
%                      7 - IQR
%                      8 - Min
%                      9 - Max
%                      10 - Range (= max - min)
%                      11 - Number of unique values (using Matlab's "unique")
%                      12 - Number of most frequent value
%                      13 - Number of 0's
%                      14 - smallest non-zero spacing
%                      15 - Continuity Index (proportion of non-zero 
%                                    pairwise distances)
%                      16 - Entropy (discrete version)
%                      17 - Bowley Skewness (robust version, based on 
%                                    quartiles and median)
%                      18 - 2nd L-statistic ratio (robust version of variance)
%                      19 - 3rd L-statistic ratio (robust version of skewness)
%                      20 - 4th L-statistic ratio (robust version of kurtosis)
if ipart == 23  
    doLog = 1;
    % Prepare data
    %nFeature = minFeature;
    mdat = vquant; % Truncate to shortest sample (min feature count)
    nFeature = size(mdat,1); % feature are quantiles
    logSuffix = '';
    if doLog == 1
        mdat = log10(mdat+0.1);
        logSuffix = '_log10';
    end    
    % Prepare param    
    savestr = [num2str(ipart) '_Q_margDistPlot_' num2str(nFeature) 'Features' logSuffix] ; % for inside plot
    nplot = 16;
    istat = 1; %<--- change

    paramstruct = struct('istat',istat,...
                         'nplot',nplot,...
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'legendcellstr',{legendcellstr}, ...
                         'mlegendcolor',mlegendcolor, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    
    % Plot marginal distribution plot
    fig = AH_figure(npc,npc,savestr);
    statstr = MargDistPlotSM(mdat,paramstruct); % data: dDimension x nSample
    % Update name of plot with stat performed
    savestr = [num2str(ipart) '_Q_margDistPlot_' statstr '_' num2str(nFeature) 'Features' logSuffix] ;
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end


%% Do PCA scatterplot
npc = 4; % number of PCs
doLog = 0;
mcolor = ones(nTeeth,1) * [0 0 1]; % initialize as blue (since scatter plot all data is black too)
mcolor(crackMask,:) = ones(sum(crackMask),1) * [1 0 0];
if ipart == 24
    %nFeature = maxFeature; % 1000 is the max, can choose smaller value to only consider earlier clusters
    mdat = vquant; % nFeatures x nSample (each column is a sample)
    nFeature = size(mdat,1);
    logSuffix = '';
    if doLog == 1
        mdat = log10(mdat+0.1);
        logSuffix = '_log10';
    end
    % Initial PCA
    paramstruct = struct('npc',npc,...
                   'iscreenwrite',1,...
                   'viout',[0 0 0 0 1]);
    outstruct = pcaSM(mdat,paramstruct);
    mpc = getfield(outstruct,'mpc'); % npc x nCluster
                        
    savestr = [num2str(ipart) '_Q_PCAScatPlot_' num2str(nFeature) 'Features' logSuffix] ;    
    paramstruct = struct('npcadiradd',4, ...
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'PC 1'; 'PC 2'; 'PC 3'; 'PC 4'}}, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1, ...
                         'isubpopkde',1); % partition data into subpopulations

    fig = AH_figure(npc,npc,savestr);
    scatplotSM(mdat,[],paramstruct);
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%% DWD analysis
if ipart == 31
    nFeature = medFeature; % 1000 is the max, can choose smaller value to only consider earlier clusters
    trainp = dataMat(1:nFeature,crackMask); % nFeatures x nSample (each column is a sample)
    trainn = dataMat(1:nFeature,~crackMask); % d x np
    mdat   = dataMat(1:nFeature,:);
    [dirvec,beta,dr] = DWD2XQ(trainp,trainn); % default penalty, handle inbalanced data
    
    savestr = [num2str(ipart) '_DWD_PCAScatPlot_' num2str(nFeature) 'Features_2PC'] ;    
    paramstruct = struct('npcadiradd',-2, ... % "-" to get ortho
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'DWD Direction'; 'OPC1'; 'OPC2'}}, ...% 'OPC2'; 'OPC3'; 'OPC3'}}, ...
                         'savestr',savestr, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'iscreenwrite',1) ;
    fig = AH_figure(npc,npc,savestr);
    scatplotSM(mdat,dirvec,paramstruct);
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
    
    savestr = [num2str(ipart) '_DWD_PCAScatPlot_' num2str(nFeature) 'Features_4PC'] ;    
    paramstruct = struct('npcadiradd',-4, ... % "-" to get ortho
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'DWD Direction'; 'OPC1'; 'OPC2'; 'OPC3'; 'OPC4'}}, ...
                         'savestr',savestr, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'iscreenwrite',1) ;
    fig = AH_figure(npc,npc,savestr);
    scatplotSM(mdat,dirvec,paramstruct);
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end


%% Do DiProPerm test on DWD class seperation
if ipart == 32   
    mlegendcolor = [1 0 0; 0 0 1]; % r, b
    savestr = [num2str(ipart) '_DiProPerm_' num2str(nFeature) 'Features'] ;    
    fig = AH_figure(2,2,savestr);
    paramstruct = struct('idir',1, ... % DWD direction vector
                         'icolor',mlegendcolor, ... % only need 2 rows of color
                         'titlecellstr',titlecellstr, ...
                         'legendcellstr',{legendcellstr}, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    [stat,pval,zscore] = DiProPermSM(trainp,trainn,paramstruct) ;
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%% DWD for quantiles
if ipart == 33
    mcolor = ones(nTeeth,1) * [0 0 1]; % initialize as blue (since scatter plot all data is black too)
    mcolor(crackMask,:) = ones(sum(crackMask),1) * [1 0 0];
    mdat = vquant; % nFeatures x nSample (each column is a sample)
    nFeature = size(mdat,1);
    trainp = mdat(1:nFeature,crackMask); % nFeatures x nSample (each column is a sample)
    trainn = mdat(1:nFeature,~crackMask); % d x np
    [dirvec,beta,dr] = DWD2XQ(trainp,trainn); % default penalty, handle inbalanced data
    
    savestr = [num2str(ipart) '_DWD_PCAScatPlot_' num2str(nFeature) 'Features_2PC'] ;    
    paramstruct = struct('npcadiradd',-2, ... % "-" to get ortho
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'DWD Direction'; 'OPC1'; 'OPC2'}}, ...% 'OPC2'; 'OPC3'; 'OPC3'}}, ...
                         'savestr',savestr, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'iscreenwrite',1) ;
    fig = AH_figure(npc,npc,savestr);
    scatplotSM(mdat,dirvec,paramstruct);
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
    
    savestr = [num2str(ipart) '_DWD_PCAScatPlot_' num2str(nFeature) 'Features_4PC'] ;    
    paramstruct = struct('npcadiradd',-4, ... % "-" to get ortho
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'DWD Direction'; 'OPC1'; 'OPC2'; 'OPC3'; 'OPC4'}}, ...
                         'savestr',savestr, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'iscreenwrite',1) ;
    fig = AH_figure(npc,npc,savestr);
    scatplotSM(mdat,dirvec,paramstruct);
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end


%% Do DiProPerm test on DWD class seperation
if ipart == 34    
    mdat = vquant; % nFeatures x nSample (each column is a sample)
    nFeature = size(mdat,1);
    trainp = mdat(1:nFeature,crackMask); % nFeatures x nSample (each column is a sample)
    trainn = mdat(1:nFeature,~crackMask); % d x np
    mlegendcolor = [1 0 0; 0 0 1]; % r, b
    savestr = [num2str(ipart) '_DiProPerm_' num2str(nFeature) 'Features'] ;    
    fig = AH_figure(2,2,savestr);
    paramstruct = struct('idir',1, ... % DWD direction vector
                         'icolor',mlegendcolor, ... % only need 2 rows of color
                         'titlecellstr',titlecellstr, ...
                         'legendcellstr',{legendcellstr}, ...
                         'isubpopkde', 1,...% partition data into subpopulations
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    [stat,pval,zscore] = DiProPermSM(trainp,trainn,paramstruct) ;
    saveas(fig, [savestr '.fig']);
    saveas(fig, [savestr '.png']);
end

%%
cd('../'); % change back to script folder