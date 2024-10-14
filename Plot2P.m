clear

% directory finding or entering code
% for now, explicitly load info
cd('/Users/sammehan/Documents/Wehr Lab/Alzheimers2P/Suite2POutput/wehr2867/10-08-24-000')


currDir = '/Users/sammehan/Documents/Wehr Lab/Alzheimers2P/Suite2POutput/wehr2867/10-08-24-000'; %maybe later this will fill with argin dir

behaviorH5 = dir('wehr*.mat');
load(behaviorH5.name)
load('Fall.mat')
fullPathH5 = fullfile(currDir, behaviorH5.name); % Need to add the ability to add a-z labels

tones = h5read('/Users/sammehan/Documents/Wehr Lab/Alzheimers2P/Suite2POutput/wehr2867/10-08-24-000/wehr2867_am_tuning_curve_20241008a.h5', '/resultsData/currentFreq');
allTones = unique(tones);

% tunedList = [14 18 24 28 30 32 35 51 58 69 72 75 84 90 96 98 104 118 135 140 161 175 187 211 213 241 259 267 301 359];

frames = info.frame;
if ~(length(frames)/2 == length(tones))
    frames = frames(1:(end-2));
end
frameIndex = 1:2:length(frames);
frames = frames(frameIndex);

for iFreq = 1:length(allTones)
    timestamps{iFreq} = frames(tones == allTones(iFreq));
end

iscellLog = logical(iscell(:, 1));
cellsToPlot = F(iscellLog, :);
neucellsToPlot = Fneu(iscellLog, :);
spikesToPlot = spks(iscellLog, :);

corrScalar = 0.7;
cellsToPlotCorr = cellsToPlot - (neucellsToPlot * corrScalar);
cellsMean = mean(cellsToPlot, 2);

for iCell = 1:sum(iscellLog)
    dFFcells(iCell, :) = (cellsToPlotCorr(iCell, :) - cellsMean(iCell))/cellsMean(iCell);
end

cmap = jet(63);

% for iCell = 1:length(tunedList)
%     currCell = tunedList(iCell);
for iCell = 1:size(dFFcells, 1)
    currCell = iCell;
    for iTone = 1:length(timestamps)
        currTimestamps = timestamps{iTone};
        if length(currTimestamps) == 63
            currTimestamps = currTimestamps(1:62);
        end
        for iTrial = 1:length(currTimestamps)
            currRange = (currTimestamps(iTrial) - 10):(currTimestamps(iTrial) + 20);
            if ~isempty(currRange(currRange <= 0))
                currRangeLog = currRange < 1;
                currRange(currRangeLog) = 1;
            end
            meanRange(iTrial, :) = dFFcells(currCell, currRange);
%             plot(spikesToPlot(currCell, currRange), 'LineWidth', 2, 'Color', cmap(iTrial,:));
%             hold on
        end
%         fig = gcf;
%         xlim([1, 31])
%         xlabel('Time (in samples, 15.49 Hz)')
% %         ylabel('dF/F')
%         xline(11, 'LineWidth', 1.5)
%         ylim([0, inf])
%         title(sprintf('ROI %s Deconvolved Spikes - %s Hz', num2str(currCell), num2str(allTones(iTone))));
%         close all

        meanTrace = mean(meanRange, 1);
        plot(meanTrace, 'k', 'LineWidth', 2)
        xlim([1, 31])
        xlabel('Time (in samples, 15.49 Hz)')
        ylabel('Mean dF/F')
        xline(11, 'LineWidth', 1.5)
        ylim([-1.5, 10])
        title(sprintf('ROI %s Mean Trace - %s Hz', num2str(currCell), num2str(allTones(iTone))));
        clear meanRange
        close all
    end
end
        
        