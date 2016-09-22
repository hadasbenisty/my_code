function [X, expLabel, combLabels, behaveDataProcessed, bahavetitles] = loadNeuronsData(datapth, files, nt, behavePath)

behaveDataProcessed = [];
bahavetitles=[];
if nargin == 3
    behavePath = '';
    bahavetitles = [];
end
dat = cell(length(files),1);
labels = cell(length(files),1);
labelNames = cell(length(files),1);
for n=1:length(files)
    labels{n} = load(fullfile(datapth, [files{n} 'sensors_dat_titles.mat']));
    labelNames{n} = getLabelNames(labels{n}.sensors_dat_titles);
    strct{n} = load(fullfile(datapth, [files{n} 'matrix.mat']));
    dat{n} = strct{n}.matrix;
    if ~isempty(behavePath)
        [behaveData{n}, bahavetitles{n}] = getBehavData(fullfile(behavePath, files{n}));
    end
end

combLabels = labelNames{1};
[X, expLabel, combLabels, behaveDataProcessed]  = getDataByCombNeurons(behavePath, files, dat, combLabels, labelNames, nt);

end


function labelNamesStr = getLabelNames(sensors_dat_titles)

for ind = 1:length(sensors_dat_titles)
    inds = strfind(sensors_dat_titles{ind}, ':');
    %     labelNames(ind) = str2double(sensors_dat_titles{ind}(inds(1)+1:inds(2)-2));
    labelNamesStr{ind} = sensors_dat_titles{ind}(inds(1)+1:inds(2)-2);
end
% [~, sortedinds] = sort(labelNames, 'ascend');
% labelNamesStr = {labelNamesStr{sortedinds}};
end



function [behaveData, bahavetitles] = getBehavData(behavePath)
s = load(fullfile(behavePath, 'points_dat_scores'));
behaveData = s.points_dat_scores;
s = load(fullfile(behavePath, 'points_dat_score_titles'));
bahavetitles = s.points_dat_score_titles;
%
% behaveData = [];
% bahavetitles = [];
% if ~isempty(behavePath)
% filesBDA = dir([behavePath '/BDA*']);
% for n = 1:length(filesBDA)
%     curr = load( fullfile(behavePath, filesBDA(n).name));
%     for m = 1:length(curr.strEvent)
%         if curr.strEvent{m}.Active
%             behavind = find(strcmp(bahavetitles, curr.strEvent{m}.Name));
%             if isempty(behavind)
%                 bahavetitles{end+1} = curr.strEvent{m}.Name;
%                 behavind = length(bahavetitles);
%             end
%             behaveData(n, behavind, m) = 1;
%         end
%     end
% end
% end
end