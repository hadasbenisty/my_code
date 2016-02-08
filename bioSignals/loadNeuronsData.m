function [X, expLabel, combLabels] = loadNeuronsData(datapth, files, nt)

dat = cell(length(files),1);
labels = cell(length(files),1);
labelNames = cell(length(files),1);
for n=1:length(files)
    labels{n} = load(fullfile(datapth, [files{n} 'sensors_dat_titles.mat']));
    labelNames{n} = getLabelNames(labels{n}.sensors_dat_titles);
    dat{n} = load(fullfile(datapth, [files{n} 'matrix.mat']));
end

combLabels = labelNames{1};

for n=2:length(files)   
    combLabels = intersect(combLabels, labelNames{n});
end

nT(1) = size(dat{1}.matrix, 2) / nt;
combinds = getCombInds(labelNames{1}, combLabels);
for T = 1:nT(1)
    X(:, :,  T) = dat{1}.matrix(combinds, 1 + nt*(T-1):nt * T);
end
expLabel = ones(nT(1), 1);
for n=2:length(files)
    nT(n) = size(dat{n}.matrix, 2) / nt;
    combinds = getCombInds(labelNames{n}, combLabels);
    for T = 1:nT(n)
        X = cat(3, X, dat{n}.matrix(combinds, 1 + nt*(T-1):nt * T));
    end
    expLabel = [expLabel; n*ones(nT(2), 1)];
end




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

function combinds = getCombInds(labelNames, combLabels)
combinds = [];
for n = 1:length(combLabels)
    if ~isempty(find(strcmp(combLabels{n}, labelNames)))
    combinds(end+1) = find(strcmp(combLabels{n}, labelNames));
    end
end
end