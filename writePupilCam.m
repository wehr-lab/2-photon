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

videoObj = VideoWriter(savename, 'MPEG-4');
videoObj.FrameRate = 16; % Should this be hard-coded or determined via file information?
open(videoObj);
writeVideo(videoObj, data);
close(videoObj);
