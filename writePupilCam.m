function [] = writePupilCam(varargin)

% Convert an *_eye.mat data structure into a MP4 Video

if ~isempty(varargin)
    datadir = convertCharsToStrings(varargin{1});
end

fullTerm = fullfile(datadir, '*_eye.mat');
matDir = dir(fullTerm);
if isempty(matDir)
    error('Could not find pupil camera data in this directory')
end
tempstr = strsplit(matDir.name, '.');
savename = fullfile(matDir.folder, tempstr{1});

pupilCam = fullfile(matDir.folder, matDir.name);

load(pupilCam);

data = data(:,:,1,(24:end-8)); %These numbers are consistent for a given requested FPS

% Uncomment this loop to determine what are the first and last frames with pupil data, and crop the video accordingly (i.e. if you change requested FPS from 24)
% for iFrame = 1:30
%     iFrameEnd = 30-iFrame; 
%     endIndex = size(data, 4) - iFrameEnd;
%     meanIntensityInitial(iFrame) = mean(data(:,:,1,iFrame),'all');
%     meanIntensityEnd(iFrame) = mean(data(:,:,1,endIndex),'all');
% end
% recStartLog = meanIntensityInitial >= 50; 
% startLog = find(recStartLog);
% recStopLog = meanIntensityEnd >= 50; 
% endLog = find(recStopLog);
% endFrame = 30-endLog(end);
% data = data(:,:,1,(startLog(1):end-endFrame));

videoObj = VideoWriter(savename, 'MPEG-4');
videoObj.FrameRate = 17.64; % Should this be hard-coded or determined via file information? 17.64 is the calculated FPS when requesting 24
open(videoObj);
writeVideo(videoObj, data);
close(videoObj);
