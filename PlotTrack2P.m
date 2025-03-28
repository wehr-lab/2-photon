function [] = PlotTrack2P(varargin)

% Pass matched_suite2p directory (/Volumes/Projects/2P5XFAD/JarascopeData/[MOUSEID]/track2p/[TRACK2P IDENTIFIER]/matched_suite2p/[SESSION ID]) 
% to plot cells tracked across sessions

if ~isempty(varargin)
    datadir = convertCharsToStrings(varargin{1});
else
    error('Gotta pass a directory to plot!')
end

datapathparts = strsplit(datadir, '/');
mouseID = datapathparts{6}; 
figdir = '/Users/sammehan/Documents/Wehr Lab/Alzheimers2P/Figs'; % where would you like to save these figures?
basedir = '/Volumes/Projects/2P5XFAD/JarascopeData/'; % full data directory path to build subsequent filepaths from
savename = fullfile(figdir, sprintf('%s-Tracked-%s.pdf', mouseID, datapathparts{8}));

matched_sessions = dir(datadir); 
matched_sessions = matched_sessions(3:end);
for iSession = 1:length(matched_sessions)
    Sessions{iSession} = matched_sessions(iSession).name;
end

for iDir = 1:length(Sessions)
    curr_session = fullfile(basedir, mouseID);
    
    tempH5 = dir(fullfile(curr_session, Sessions{iDir}, '*.h5'));
    H5Paths{iDir} = fullfile(basedir, mouseID, Sessions{iDir}, tempH5.name);
    tempMat = dir(fullfile(curr_session, Sessions{iDir}, '*.mat'));
    MatPaths{iDir} = fullfile(basedir, mouseID, Sessions{iDir}, tempMat.name);
    
    tempFall = dir(fullfile(datadir, Sessions{iDir}, '/suite2p/plane0/Fall.mat'));
    if isempty(tempFall)
        F = readNPY(fullfile(datadir, Sessions{iDir}, '/suite2p/plane0/F.npy'));
        Fneu = readNPY(fullfile(datadir, Sessions{iDir}, '/suite2p/plane0/Fneu.npy'));
        iscell = readNPY(fullfile(datadir, Sessions{iDir}, '/suite2p/plane0/iscell.npy'));
        spks = readNPY(fullfile(datadir, Sessions{iDir}, '/suite2p/plane0/spks.npy'));
%         stat = readNPY(fullfile(datadir, Sessions{iDir}, '/suite2p/plane0/stat.npy'));
%         ops = readNPY(fullfile(datadir, Sessions{iDir}, '/suite2p/plane0/ops.npy'));
%         save(fullfile(tempFall, 'Fall.mat'), 'F', 'Fneu', 'iscell', 'ops', 'spks', 'stat');
%         clear F Fneu iscell ops spks stat
    end
%     FallPaths{iDir} = fullfile(tempFall.folder, 'Fall.mat');
    tones = h5read(H5Paths{iDir}, '/resultsData/currentFreq');
    intensities = h5read(H5Paths{iDir}, '/resultsData/currentIntensity');
    allTones = unique(tones);
    allInts = unique(intensities); allInts = flip(allInts);
    
    load(MatPaths{iDir})
    frames = info.frame;
    if rem(length(info.frame), length(tones)) == 2
    elseif ~(length(frames)/2 == length(tones))
        frames = frames(1:(end-2));
    end
    frameIndex = 1:2:length(frames);
    frames = frames(frameIndex);
    
    nCond = 0;
    for iFreq = 1:length(allTones)
        for iInt = 1:length(allInts)
            nCond = nCond + 1;
            tempTimestamps = frames(tones == allTones(iFreq));
            tempTimestampsInt = frames(intensities == allInts(iInt));
            timestamps{iFreq, iInt} = tempTimestamps(ismember(tempTimestamps, tempTimestampsInt));
            nReps(nCond) = length(timestamps{iFreq, iInt});
        end
    end
    minReps = min(nReps);
    allTimestamps{iDir} = timestamps;
    allMinReps(iDir) = minReps;
    clear timestamps 
    
    if exist('iscell') == 1
        iscellList = iscell;
        clear iscell
    end
    iscellLog = logical(iscellList(:, 1)); iscellThresh = iscellList(:, 2);
    % Uncomment and enter a threshold value (0-1) to use Suite2P's likelihood value to select good cells
    % S2Pthresh = 0.95; % Suite2P likelihood threshold to use
    % iscellLog = iscellThresh >= S2Pthresh; 
    cellsToPlot = F(iscellLog, :);
    neucellsToPlot = Fneu(iscellLog, :);

    corrScalar = 0.7;
    cellsToPlotCorr{iDir} = cellsToPlot - (neucellsToPlot * corrScalar);
end    

for currCell = 1:size(cellsToPlotCorr{1}, 1)
    for iDir = 1:length(Sessions)
        timestamps = allTimestamps{iDir};
        for iTone = 1:size(timestamps, 1)
            for iInt = 1:size(timestamps, 2)
                
                currTimestamps = timestamps{iTone, iInt};
                if length(currTimestamps) > allMinReps(iDir) || iDir == 1
                    currTimestamps = currTimestamps(1:allMinReps(iDir));
                end
                for iTrial = 1:length(currTimestamps)
                    currRange = (currTimestamps(iTrial) - 10):(currTimestamps(iTrial) + 20);
                    normRange = (currTimestamps(iTrial) - 11):(currTimestamps(iTrial) -1);
                    if ~isempty(currRange(currRange <= 0))
                        currRangeLog = currRange < 1;
                        currRange(currRangeLog) = 1;
                    end
                    if ~isempty(currRange(currRange > frames(end)))
                        currRangeLog = currRange > size(cellsToPlotCorr{iDir}, 2);
                        currRange(currRangeLog) = size(cellsToPlotCorr{iDir}, 2);
                    end
                    if sum(normRange <= 0) == length(normRange) || sum(normRange <= 0) >= 1
                        normRange = (currTimestamps(iTrial) + 21):(currTimestamps(iTrial) + 30);
                    end
                    currTrace = (cellsToPlotCorr{iDir}(currCell, currRange) - mean(cellsToPlotCorr{iDir}(currCell, normRange)))/mean(cellsToPlotCorr{iDir}(currCell, normRange));
                    meanRange(iTrial, :) = currTrace;
                end
                meanRanges{iTone, iInt} = meanRange';
            end
        end

        if iDir == 1
            subplot1(2,6, 'Min', [0.05, 0.05], 'Gap', [0.01, 0.01]);
            fig = gcf; orient(fig, 'landscape');
            axes(fig, 'Position', [0.05, 0.05, 0.9, 0.9])
            title(sprintf('ROI %s Tuning Curve', num2str(currCell)), 'Position', [0.5, 1.02]);
            text(0.5, -0.04, 'Time (in samples, 15.49 Hz)', 'HorizontalAlignment', 'center');
            ylabel('dF/F')
            axis off
            gcf;
        else
            figure
            fig = gcf; orient(fig, 'landscape');
            axes(fig, 'Position', [0.05, 0.05, 0.9, 0.9])
            axis off
            gcf;
        end

        which_fig = 0;
        for iInt = 1:size(timestamps, 2)
            for iTone = 1:size(timestamps, 1)
                which_fig = which_fig + 1;
                if sum(isnan(meanRanges{iTone, iInt}), 'all') ~= 0
                    meanRanges{iTone, iInt} = rmmissing(meanRanges{iTone, iInt}, 2);
                end
                meanTrace = mean(meanRanges{iTone, iInt}, 2);
                if iDir == 1
                    subplot1(which_fig);
                end
                plot(meanRanges{iTone, iInt}, 'r', 'LineWidth', 1); hold on; plot(meanTrace, 'k', 'LineWidth', 2); xlim([1, 31]); xline(11, 'k', 'LineWidth', 1.5); ylim([-1.5, 10]);
                if which_fig == 1 && iDir == 1
                    ylabel('dF/F - 70 dbSPL');
                elseif which_fig == 7 && iDir == 1
                    ylabel('dF/F - 50 dbSPL');
                elseif iDir == 2
                    xlabel('Time (in samples, 15.49 Hz)', 'HorizontalAlignment', 'center');
                    ylabel('dF/F')
                    title(sprintf('ROI %s WN Response', num2str(currCell)));
                else
                end
                if which_fig == 1 && iDir == 1
                    xlabel('2000 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 2 && iDir == 1
                    xlabel('3482 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 3 && iDir == 1
                    xlabel('6063 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 4 && iDir == 1
                    xlabel('10556 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 5 && iDir == 1
                    xlabel('18379 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 6 && iDir == 1
                    xlabel('32000 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                end
                clear meanTrace
            end
        end
        %print(savename, '-dpsc2', '-append', '-bestfit');
        if currCell == 1
            exportgraphics(gcf, savename);
        else
            exportgraphics(gcf, savename, 'Append', true);
        end
        sprintf('On Cell %d / %d \n', currCell, size(cellsToPlotCorr{iDir}, 1))
        clear meanRange meanRanges normRange currTrace
        close all
    end
end