function fig = Plot2PCurves(plotDir, iCell)

plotDeconvolved = 0; %logical for plotting deconvolved signal
currDir = plotDir; % enter directory containing .sbx and .mat recording files
% tunedList = [14 18 24 28 30 32 35 51 58 69 72 75 84 90 96 98 104 118 135 140 161 175 187 211 213 241 259 267 301 359];

cd(currDir)

matH5 = dir('wehr*.mat');
load(matH5.name)
FallDir = dir('Fall.mat');
if isempty(FallDir)
    currDirContents = dir('suite2p');
    if ~isempty(currDirContest)
        s2pDirContents = dir('~/suite2p/plane0');
        if ~isempty(s2pDirContents)
            plane0DirContents = dir('~/suite2p/plane0/Fall.mat');
            if ~isempty(plane0DirContents)
                load('Fall.mat')
                iscellList = load('Fall.mat', 'iscell');
                iscellList = iscellList.iscell(:, 1);
                clear iscell
            else
                error('Could not locate F/all statistics in this directory!')
            end
        else
            error('Could not locate F/all statistics in this directory!')
        end
    else
        error('Could not locate F/all statistics in this directory!')
    end
else
    load('Fall.mat')
    iscellList = load('Fall.mat', 'iscell');
    iscellList = iscellList.iscell(:, 1);
    clear iscell
end

behaviorH5 = dir(fullfile(plotDir, 'wehr*.h5'));
if isempty(behaviorH5)
    currDirContents = dir('suite2p');
    if ~isempty(currDirContents)
        cd('suite2p')
        s2pDirContents = dir('plane0');
        if ~isempty(s2pDirContents)
            cd('plane0')
            plane0DirContents = dir('wehr*.h5');
            if ~isempty(plane0DirContents)
                fullpathH5 = fullfile(pwd, plane0DirContents.name);
            else
                error('Could not locate H5 Behavior file, make sure you transferred it over from the "behavior" repository!')
            end
        end
    end
else
    fullpathH5 = fullfile(currDir, behaviorH5.name);
end

tones = h5read(fullpathH5, '/resultsData/currentFreq');
allTones = unique(tones);

frames = info.frame;
if ~(length(frames)/2 == length(tones))
    frames = frames(1:(end-2));
end
frameIndex = 1:2:length(frames);
frames = frames(frameIndex);

for iFreq = 1:length(allTones)
    timestamps{iFreq} = frames(tones == allTones(iFreq));
end

iscellLog = logical(iscellList(:, 1));
cellsToPlot = F(iscellLog, :);
neucellsToPlot = Fneu(iscellLog, :);
spikesToPlot = spks(iscellLog, :);

corrScalar = 0.7;
cellsToPlotCorr = cellsToPlot - (neucellsToPlot * corrScalar);
cellsMean = mean(cellsToPlot, 2);

for i = 1:length(timestamps)
    timestampCount(i) = length(timestamps{1, i});
end
minTimestamps = min(timestampCount);
cmap = jet(max(timestampCount));

currCell = iCell;
for iTone = 1:length(timestamps)
    currTimestamps = timestamps{iTone};
    if length(currTimestamps) ~= minTimestamps
        currTimestamps = currTimestamps(1:minTimestamps);
    end
    for iTrial = 1:length(currTimestamps)
        currRange = (currTimestamps(iTrial) - 10):(currTimestamps(iTrial) + 20);
        normRange = (currTimestamps(iTrial) - 11):(currTimestamps(iTrial) - 1);
        if ~isempty(currRange(currRange <= 0))
            currRangeLog = currRange < 1;
            currRange(currRangeLog) = 1;
        end
        if ~isempty(normRange(normRange <= 0))
            normRange = normRange(normRange > 0);
        end
        currTrace = (cellsToPlotCorr(currCell, currRange) - mean(cellsToPlotCorr(currCell, normRange)))/mean(cellsToPlotCorr(currCell, normRange));
        meanRange(iTrial, :) = currTrace;
        if plotDeconvolved == 1
            if iTrial == 1
                figDecon = figure;
            end
            plot(spikesToPlot(currCell, currRange), 'LineWidth', 2, 'Color', cmap(iTrial,:));
            hold on
        end
    end
    if plotDeconvolved == 1
        xlim([1, 31])
        xlabel('Time (in samples, 15.49 Hz)')
        ylabel('Deconvolved Spike Signal')
        xline(11, 'LineWidth', 1.5)
        ylim([0, inf])
        title(sprintf('ROI %s Deconvolved Spikes - %s Hz', num2str(currCell), num2str(allTones(iTone))));
        hold off
    end
    
    fig = figure;
    meanTrace = mean(meanRange, 1);
    hold on
    plot(meanRange', 'k', 'LineWidth', 1)
    plot(meanTrace, 'r', 'LineWidth', 2)
    xlim([1, 31])
    xlabel('Time (in samples, 15.49 Hz)')
    ylabel('Mean dF/F')
    xline(11, 'LineWidth', 1.5)
    ylim([-1.5, 10])
    title(sprintf('ROI %s Mean Trace - %s Hz', num2str(currCell), num2str(allTones(iTone))));
    clear meanRange
end