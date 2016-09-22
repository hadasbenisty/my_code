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
addpath('newcodeuri');
%%%
% Parameters
%%%
% directory where TPA and BDA files located (Change it if you need)
pth = '..\..\..\datasets\biomed\D30new';
pth = '..\..\..\datasets\biomed\D16';

dirs = dir(pth);
for n=11:length(dirs)
    analysisDir         = fullfile('..\..\..\datasets\biomed\D16\', dirs(n).name);
    
    % frame rate ratio between Two Photon imaging and Behavior
    frameRateRatio      = 6;%6; % (do not touch)
    % Common names to be converted to Ids
    eventNameList      = {'Lift','Grab','Sup','Atmouth','Chew','Sniff','Handopen','Botharm','Tone','Table','Failure','Success', 'Trying'}; % (do not touch)
    
    %%%
    % Load Two Photon ROI data
    %%%
    fileNames        = dir(fullfile(analysisDir,'TPA_*.mat'));
    fileNum          = length(fileNames);
    if fileNum < 1,
        error('TwoPhoton : Can not find data files in the directory %s. Check file or directory names.',analysisDir);
    end;
    [fileNamesRoi{1:fileNum,1}] = deal(fileNames.name);
    fileNumRoi                  = fileNum;
    
    allTrialRois                    = cell(fileNumRoi,1);
    for trialInd = 1:fileNumRoi,
        fileToLoad                 = fullfile(analysisDir,fileNamesRoi{trialInd});
        usrData                    = load(fileToLoad);
        allTrialRois{trialInd}     = usrData.strROI;
    end
    
    % match new format to old format, the deltaF over F is saved in Data(:,2)
    % instead of procROI
    if ~isfield(allTrialRois{trialInd}{1}, 'Data') && ~isprop(allTrialRois{trialInd}{1}, 'Data')
        if isfield(allTrialRois{trialInd}{1}, 'procROI')
            for trialInd = 1:fileNumRoi,
                for m = 1:length(allTrialRois{trialInd})
                    allTrialRois{trialInd}{m}.Data(:,2) = allTrialRois{trialInd}{m}.procROI;
                end
            end
        else
            error('Unclear on which field the data is saved');
        end
    end
    %%%
    % Load Behavioral Event data
    %%%
    fileNames        = dir(fullfile(analysisDir,'BDA_*.mat'));
    fileNum          = length(fileNames);
    if fileNum < 1,
        error('Behavior : Can not find data files in the directory %s. Check file or directory names.',analysisDir);
    end;
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
    roisPerTrialNum   = length(allTrialRois{1});    
    dffDataArray     = zeros(roisPerTrialNum,length(allTrialRois{1}{1}.Data(:,2)),length(allTrialRois));
    
    for trial_i = 1:min(fileNumRoi,fileNumEvent)
        %%%
        % Extract specific trial data for Two Photon ROI
        %%%
        % check
        if trial_i < 1 || trial_i > fileNumRoi,
            error('Requested trial should be in range [1:%d]',fileNumRoi)
        end
        if trial_i < 1 || trial_i > fileNumEvent,
            error('Requested trial should be in range [1:%d]',fileNumEvent)
        end
        
        % extract ROI dF/F data
        
        if roisPerTrialNum < 1,
            error('Could not find ROI data for trial %d',trial_i)
        end
        dffData          = allTrialRois{trial_i}{1}.Data(:,2); % actual data
        [framNum,~]      = size(dffData);
        
        roiNames         = cell(1,roisPerTrialNum);
        % collect
        
        for m = 1:roisPerTrialNum,
            if isempty(allTrialRois{trial_i}{m}.Data),
                error('Two Photon data is not complete. Can not find dF/F results.')
            end
            dffDataArray(m,:, trial_i) = allTrialRois{trial_i}{m}.Data(:,2);
            roiLoc{m} = allTrialRois{trial_i}{m}.xyInd;
            roiNames{m}       = allTrialRois{trial_i}{m}.Name;
        end
    end
    
    %%%
    % Extract specific trial data for Behavioral Events
    %%%
    % extract Event Time data
    eventsPerTrialNum = 0;
    for trial_i = 1:length(allTrialRois)
        eventsPerTrialNum   = max(eventsPerTrialNum, length(allTrialEvents{trial_i}));
    end
    if eventsPerTrialNum < 1,
        warning('Could not find Event data for trial %d',trial_i);
    end
    timeData         = allTrialEvents{trial_i}{1}.TimeInd; % actual data
    eventDataArray   = zeros(framNum,eventsPerTrialNum, length(allTrialRois));
    eventNames       = cell(1,eventsPerTrialNum);
    
    % collect
    for trial_i = 1:length(allTrialRois)
        for m = 1:length(allTrialEvents{trial_i}),
            timeInd     = allTrialEvents{trial_i}{m}.tInd;
            timeInd     = round(timeInd./frameRateRatio); % transfers to time of the two photon
            timeInd     = max(1,min(framNum,timeInd));
            % assign to vector
            eventDataArray(timeInd(1):timeInd(2),m, trial_i) = allTrialEvents{trial_i}{m}.Id;
            eventNames{m} = sprintf('%s - %d',allTrialEvents{trial_i}{m}.Name, allTrialEvents{trial_i}{m}.SeqNum);
        end
    end

%%%
% Show
%%%
timeImage       = 1:framNum;
figure,
subplot(221)
plot(timeImage,dffDataArray(:,:,trial_i)),legend(roiNames)
xlabel('Time [Frame]'), ylabel('dF/F'), title(sprintf('Two Photon data for trial %d',trial_i))

subplot(222)
plot(timeImage,eventDataArray(:,:,trial_i)),legend(eventNames)
xlabel('Time [Frame]'), ylabel('Events'), title(sprintf('Event duration data for trial %d',trial_i))

subplot(223);imagesc(dffDataArray(:,:,trial_i))
subplot(224);imagesc(eventDataArray(:,:,trial_i).')


save(fullfile(analysisDir, 'sampsAndBehave'), 'dffDataArray', 'eventDataArray', 'roiNames', 'eventNameList', 'roiLoc');
end
return
