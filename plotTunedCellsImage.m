function [fig, fig2] = plotTunedCellsImage(plotDir)

cd(plotDir)
matOps = dir('wehr*.mat');
load(matOps.name);
load('Fall.mat')

iscellList = load('Fall.mat', 'iscell');
iscellList = iscellList.iscell;
clear iscell

iscellLog = logical(iscellList(:, 1));
cellsToPlot = F(iscellLog, :);
neucellsToPlot = Fneu(iscellLog, :);
spikesToPlot = spks(iscellLog, :);

h5Ops = dir('wehr*.h5');

h5path = fullfile(plotDir, h5Ops.name);

tones = h5read(h5path, '/resultsData/currentFreq');
allTones = unique(tones);

goodStats = stat(iscellLog);
%tunedList = [14 18 24 28 30 32 35 51 58 69 72 75 84 90 96 98 104 118 135 140 161 175 187 211 213 241 259 267 301 359];
tunedList = [3 13 18 32 53 59 71 73 114 120 125 289 306 307];


for iROI = 1:length(tunedList)
    if iROI == 1
        currStats = goodStats(tunedList);
        fig = figure;
    end
    currStats{1, iROI}.ypix = abs(currStats{1, iROI}.ypix - 512);
    scatter(currStats{1, iROI}.xpix, currStats{1, iROI}.ypix, 8, '.', 'LineWidth', 8, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    hold on
end
xlim([0 512]);
ylim([0 512]);
hold off
fig2 = figure;
imagesc(ops.meanImg);