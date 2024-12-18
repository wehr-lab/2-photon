function [] = Plot2P(varargin)

    % You can either pass the full path to your 2P datadir (on the NAS), or you can run as a script by entering your plot directory below
    
    if ~isempty(varargin)
        datadir = varargin{1};
    else
        datadir = '/Volumes/Projects/2P5XFAD/JarascopeData/wehr3133/12-12-24-000'; % enter directory to plot (if none explicitly passed)
    end
    cd(datadir)

    figdir = '/Users/sammehan/Documents/Wehr Lab/Alzheimers2P/Figs';
    temp1 = strsplit(datadir, '/');
    savename = fullfile(figdir, strcat(temp1{end-1}, '-', temp1{end}, '-TC.ps'));

    behaviorMAT = dir('wehr*.mat');
    load(behaviorMAT.name)
    FallPath = fullfile(datadir, '/suite2p/plane0/Fall.mat');
    load(FallPath)
    iscellList = load(FallPath, 'iscell');
    iscellList = iscellList.iscell;
    clear iscell
    fullPathMAT = fullfile(datadir, behaviorMAT.name); % Need to add the ability to add a-z labels, or we just presort each behavior h5 with the correct Ephys dir

    behaviorH5 = dir('wehr*.h5');
    fullPathH5 = fullfile(datadir, behaviorH5.name);

    tones = h5read(fullPathH5, '/resultsData/currentFreq');
    intensities = h5read(fullPathH5, '/resultsData/currentIntensity');
    times = h5read(fullPathH5, '/events/eventTime');
    codes = h5read(fullPathH5, '/events/eventCode');
    allTones = unique(tones);
    allInts = unique(intensities); allInts = flip(allInts);
    % tones = tones(837:end);
    % intensities = intensities(837:end); % +137; ignore this, manual coding to overcome TasKontrol appending data sessions

    frames = info.frame;
    if ~(length(frames)/2 == length(tones))
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

    iscellLog = logical(iscellList(:, 1));
    cellsToPlot = F(iscellLog, :);
    neucellsToPlot = Fneu(iscellLog, :);
    spikesToPlot = spks(iscellLog, :);

    corrScalar = 0.7;
    cellsToPlotCorr = cellsToPlot - (neucellsToPlot * corrScalar);
    cellsMean = mean(cellsToPlot, 2);

    for iCell = 1:size(cellsToPlotCorr, 1)
        currCell = iCell;
        for iTone = 1:length(timestamps)
            for iInt = 1:length(allInts)
                currTimestamps = timestamps{iTone, iInt};
                if length(currTimestamps) > minReps
                    currTimestamps = currTimestamps(1:minReps);
                end
                for iTrial = 1:length(currTimestamps)
                    currRange = (currTimestamps(iTrial) - 10):(currTimestamps(iTrial) + 20);
                    normRange = (currTimestamps(iTrial) - 11):(currTimestamps(iTrial) -1);
                    if ~isempty(currRange(currRange <= 0))
                        currRangeLog = currRange < 1;
                        currRange(currRangeLog) = 1;
                    end
                    if ~isempty(currRange(currRange > frames(end)))
                        currRangeLog = currRange > size(cellsToPlotCorr, 2);
                        currRange(currRangeLog) = size(cellsToPlotCorr, 2);
                    end
                    if sum(normRange <= 0) == length(normRange)
                        normRange = (currTimestamps(iTrial) + 21):(currTimestamps(iTrial) + 30);
                    end
                    currTrace = (cellsToPlotCorr(currCell, currRange) - mean(cellsToPlotCorr(currCell, normRange)))/mean(cellsToPlotCorr(currCell, normRange));
                    meanRange(iTrial, :) = currTrace;
                end
                meanRanges{iTone, iInt} = meanRange';
            end
        end

        subplot1(2,6, 'Min', [0.05, 0.05], 'Gap', [0.01, 0.01]);
        fig = gcf; orient(fig, 'landscape');
        axes(fig, 'Position', [0.05, 0.05, 0.9, 0.9])
        title(sprintf('ROI %s Tuning Curve', num2str(currCell)), 'Position', [0.5, 1.02]); 
        text(0.5, -0.04, 'Time (in samples, 15.49 Hz)', 'HorizontalAlignment', 'center');
        axis off
        gcf;

        which_fig = 0;
        for iInt = 1:length(allInts)
            for iTone = 1:length(timestamps)
                which_fig = which_fig + 1;
                if sum(isnan(meanRanges{iTone, iInt}), 'all') ~= 0
                    meanRanges{iTone, iInt} = rmmissing(meanRanges{iTone, iInt}, 2);
                end
                meanTrace = mean(meanRanges{iTone, iInt}, 2);

                subplot1(which_fig); plot(meanRanges{iTone, iInt}, 'r', 'LineWidth', 1); plot(meanTrace, 'k', 'LineWidth', 2); xlim([1, 31]); xline(11, 'LineWidth', 1.5); ylim([-1.5, 10]); 
                if which_fig == 1
                    ylabel('dF/F - 70 dbSPL');
                elseif which_fig == 7
                    ylabel('dF/F - 50 dbSPL');
                end
                if which_fig == 1
                    xlabel('2000 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 2
                    xlabel('3482 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 3
                    xlabel('6063 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 4
                    xlabel('10556 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 5
                    xlabel('18379 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                elseif which_fig == 6
                    xlabel('32000 Hz', 'Position', [11, 10.65], 'HorizontalAlignment', 'center');
                end
                clear meanTrace
            end    
        end
        print(savename, '-dpsc2', '-append', '-bestfit');
        sprintf('On Cell %d / %d \n', currCell, size(cellsToPlotCorr, 1))
        clear meanRange meanRanges normRange currTrace
        close all
    end

