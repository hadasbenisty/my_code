% TPA_ShowEventsRois
% Script to show how to use ROI and Event data structures.
% Select analysis directory that contains TPA_*.mat and BDA_*.mat files, load them and display.
% Inputs:
%       directory to be analysed
% Outputs:
%        view of the data

% Event structure description:
% strEvent.Type        = 1;             % ROI_TYPES.RECT 
% strEvent.Active      = true;          % designates if this pointer structure is in use
% strEvent.NameShow    = false;         % manage show name
% strEvent.zInd        = 1;             % location in Z stack
% strEvent.tInd        = 1;             % location in T stack
% strEvent.Position    = pos;           % rect position
% strEvent.xyInd       = xy;            % shape in xy plane
% strEvent.Name        = Name from the list and animal name;
% strEvent.TimeInd     = [start stop];  % time/frame indices
% strEvent.SeqNum      = 1;             % designates Event sequence number in one trial
% strEvent.Color       = [0 1 0];       % designates Event color


%-----------------------------
% Ver	Date	 Who	Descr
%-----------------------------
% 19.09 14.10.14 UD     Adding more events
% 19.08 11.10.14 UD     Adding enumeration to event names and sequence nummber
% 18.12 11.07.14 UD     created
%-----------------------------

%%%
% Parameters
%%%
% directory where TPA and BDA files located (Change it if you need)
% analysisDir         = 'C:\Users\mayahar\Documents\maya-research\Questionnaire\Data\D8\8_14_14\';
%analysisDir         = 'C:\LabUsers\Uri\Data\Janelia\Analysis\m2\4_4_14';
function create_control_cno_data_for_quest

controlDir     = '..\..\..\datasets\biomed\8_6_14_1-20_control';
cnoDir         = '..\..\..\datasets\biomed\8_6_14_21-60_cno';
% controlDir     = 'C:\Users\mayahar\Documents\maya-research\Questionnaire\Data\D10\8_6_14 new\8_6_14_control\';
% cnoDir         = 'C:\Users\mayahar\Documents\maya-research\Questionnaire\Data\D10\8_6_14 new\8_6_14_cno\';

% controlDir         = 'C:\Users\galga\Dropbox\research\studies\Adavanced topics in Signal Processing\D10\8_6_14_1-20_control\';
% cnoDir         = 'C:\Users\galga\Dropbox\research\studies\Adavanced topics in Signal Processing\D10\8_6_14_21-60_cno\';
% special trial (repetition) to show (Change it if you need)
trialIndShow        = 1;


%%%
% Load Two Photon ROI data - control
%%%
fileNames        = dir(fullfile(controlDir,'TPA_*.mat'));
fileNum          = length(fileNames);
if fileNum < 1,
    error('TwoPhoton : Can not find data files in the directory %s. Check file or directory names.',controlDir);
end;
[fileNamesRoi{1:fileNum,1}] = deal(fileNames.name);
controlFileNumRoi                  = fileNum;

controlTrialRois               = cell(controlFileNumRoi,1);
for trialInd = 1:controlFileNumRoi,
    fileToLoad                 = fullfile(controlDir,fileNamesRoi{trialInd});
    usrData                    = load(fileToLoad);
    controlTrialRois{trialInd} = usrData.strROI;
end

%%%
% Load Two Photon ROI data - cno
%%%
fileNames        = dir(fullfile(cnoDir,'TPA_*.mat'));
fileNum          = length(fileNames);
if fileNum < 1,
    error('TwoPhoton : Can not find data files in the directory %s. Check file or directory names.',cnoDir);
end;
[fileNamesRoi{1:fileNum,1}] = deal(fileNames.name);
cnoFileNumRoi                  = fileNum;

cnoTrialRois                   = cell(cnoFileNumRoi,1);
for trialInd = 1:cnoFileNumRoi,
    fileToLoad                 = fullfile(cnoDir,fileNamesRoi{trialInd});
    usrData                    = load(fileToLoad);
    cnoTrialRois{trialInd}     = usrData.strROI;
end

%%%%%
controlRoisPerTrial   = length(controlTrialRois{1});
controlRoiNames         = cell(1,controlRoisPerTrial);
for m = 1:controlRoisPerTrial,
    if isempty(controlTrialRois{1}{m}.procROI),
        error('Two Photon data is not complete. Can not find dF/F results.')
    end
    controlRoiNames{m}   = controlTrialRois{1}{m}.Name;
end

cnoRoisPerTrial   = length(cnoTrialRois{1});
cnoRoiNames         = cell(1,cnoRoisPerTrial);
for m = 1:cnoRoisPerTrial,
    if isempty(cnoTrialRois{1}{m}.procROI),
        error('Two Photon data is not complete. Can not find dF/F results.')
    end
    cnoRoiNames{m}       = cnoTrialRois{1}{m}.Name;
end

[roiNames,controlRoiInds,cnoRoiInds]= intersect(controlRoiNames,cnoRoiNames);
controlRoisPerTrial   = length(roiNames);
cnoRoisPerTrial       = length(roiNames);

matrix_control = get_data(controlDir, controlTrialRois, controlRoiInds, controlFileNumRoi, roiNames);
matrix_cno     = get_data(cnoDir,     cnoTrialRois,     cnoRoiInds,     cnoFileNumRoi,     roiNames);


function matrix = get_data(analysisDir, allTrialRois, roiInds, fileNumRoi, roiNames)
% frame rate ratio between Two Photon imaging and Behavior
frameRateRatio      = 6*3;%6; % (do not touch)
% Common names to be converted to Ids
eventNameList      = {'Lift','Grab','Sup','Atmouth','Chew','Sniff','Handopen','Botharm','Tone','Table','Miss'}; % (do not touch)
nEvents = length(eventNameList);

%%%
% Load Behavioral Event data
%%%
fileNames        = dir(fullfile(analysisDir,'BDA_*.mat'));
fileNum          = length(fileNames);
if fileNum < 1,
    warning('Behavior : Can not find event files in the directory %s. Check file or directory names.',analysisDir);
    eventsLabelExist = false;
    num_trials = fileNumRoi;
else
    [fileNamesEvent{1:fileNum,1}] = deal(fileNames.name);
    fileNumEvent                  = fileNum;
    
    allTrialEvents                = cell(fileNumEvent,1);
    for trialInd = 1:fileNumEvent,
        fileToLoad                 = fullfile(analysisDir,fileNamesEvent{trialInd});
        usrData                    = load(fileToLoad);
        allTrialEvents{trialInd}   = usrData.strEvent;
    end
    
    %%%
    % Convert Event names to Ids
    %%%
    % add field Id
    eventNameList                  = lower(eventNameList);
    for trialInd = 1:fileNumEvent,
        eventNum                   = length(allTrialEvents{trialInd});
        for m = 1:eventNum,
            eName                      = lower(allTrialEvents{trialInd}{m}.Name);
            eBool                      = strncmp(eName, eventNameList,3);
            eInd                       = find(eBool);
            if isempty(eInd),
                fprintf('W : Trial %d, Event %d has undefined name %s. Set Id = 0.\n',trialInd,m,allTrialEvents{trialInd}{m}.Name)
                eInd = 0;
            end
            allTrialEvents{trialInd}{m}.Id = eInd;
        end
    end
    num_trials = min(fileNumRoi,fileNumEvent);
    eventsLabelExist = true;
end

roisPerTrial   = length(roiInds);

matrix                   = []; % data structure - rows are sensors (ROIs), cols are time
if eventsLabelExist
    points_dat.scores        = []; % nTimes*nEvents binary (event is on or off)
    points_dat.score_titles  = eventNameList;
end
sensors_dat.titles       = roiNames;

for trialIndShow = 1:num_trials
    %%%
    % Extract specific trial data for Two Photon ROI
    %%%
    % check
    if trialIndShow < 1 || trialIndShow > fileNumRoi,
        error('Requested trial should be in range [1:%d]',fileNumRoi)
    end
    if eventsLabelExist && (trialIndShow < 1 || trialIndShow > fileNumEvent)
        error('Requested trial should be in range [1:%d]',fileNumEvent)
    end
    
    % extract ROI dF/F data
    
    allTrialRois{trialIndShow} = allTrialRois{trialIndShow}(roiInds);
    roisPerTrialNum   = length(allTrialRois{trialIndShow});
    if roisPerTrialNum < 1,
        error('Could not find ROI data for trial %d',trialIndShow)
    elseif roisPerTrialNum ~= roisPerTrial
        error('trial %d has more ROIs', trialIndShow)
    end
    dffData          = allTrialRois{trialIndShow}{1}.procROI; % actual data
    [framNum,~]      = size(dffData);
    dffDataArray     = repmat(dffData,1,roisPerTrialNum);
    trialRoiNames         = cell(1,roisPerTrialNum);
    % collect
    for m = 1:roisPerTrialNum,
        if isempty(allTrialRois{trialIndShow}{m}.procROI),
            error('Two Photon data is not complete. Can not find dF/F results.')
        end
        dffDataArray(:,m) = allTrialRois{trialIndShow}{m}.procROI;
        trialRoiNames{m}  = allTrialRois{trialIndShow}{m}.Name;
    end
    
    if ~all(strcmp(trialRoiNames,roiNames))
        error('Different ROI names in trial %d',trialIndShow)
    end
    
    %%%
    % Extract specific trial data for Behavioral Events
    %%%
    % extract Event Time data
    if eventsLabelExist
        eventsPerTrialNum   = length(allTrialEvents{trialIndShow});
        if eventsPerTrialNum < 1,
            warning('Could not find Event data for trial %d',trialIndShow);
        end
        timeData         = allTrialEvents{trialIndShow}{1}.TimeInd; % actual data
        eventDataArray   = zeros(framNum,eventsPerTrialNum);
        eventNames       = cell(1,eventsPerTrialNum);
        
        trialEvents = zeros(nEvents,framNum);
        % collect
        for m = 1:eventsPerTrialNum,
            timeInd     = allTrialEvents{trialIndShow}{m}.TimeInd;
            timeInd     = round(timeInd./frameRateRatio); % transfers to time of the two photon
            timeInd     = max(1,min(framNum,timeInd));
            % assign to vector
            eventDataArray(timeInd(1):timeInd(2),m) = 1;
            eventNames{m} = sprintf('%s - %d',allTrialEvents{trialIndShow}{m}.Name, allTrialEvents{trialIndShow}{m}.SeqNum);
            trialEvents(allTrialEvents{trialIndShow}{m}.Id,timeInd(1):timeInd(2)) = 1;
        end
    end
    
    %%%
    matrix = [matrix , dffDataArray(2:end,:)'];
    if eventsLabelExist
        points_dat.scores  = [points_dat.scores  ,trialEvents(:,2:end)];
    end
    
end

figure;imagesc(matrix);
sensors_dat_titles = sensors_dat.titles;
save([analysisDir 'matrix.mat'],'matrix')
save([analysisDir 'sensors_dat_titles.mat'],'sensors_dat_titles')
fileID = fopen([analysisDir 'sensors_dat_titles.txt'],'w');
formatSpec = '%s\n';
[nrows,ncols] = size(sensors_dat_titles);
for row = 1:nrows
    fprintf(fileID,formatSpec,sensors_dat_titles{row,:});
end
fclose(fileID);

if eventsLabelExist
    figure;imagesc(points_dat.scores);
    points_dat_scores = points_dat.scores;
    points_dat_score_titles = points_dat.score_titles;
    save([analysisDir 'points_dat_scores.mat'],'points_dat_scores')
    save([analysisDir 'points_dat_score_titles.mat'],'points_dat_score_titles')
    fileID = fopen([analysisDir 'points_dat_score_titles.txt'],'w');
    formatSpec = '%s\n';
    [nrows,ncols] = size(points_dat_score_titles);
    for row = 1:nrows
        fprintf(fileID,formatSpec,points_dat_score_titles{row,:});
    end
    fclose(fileID);
end

% save -ascii 'quest_matrix.txt' matrix
return