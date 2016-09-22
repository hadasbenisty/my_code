function behaveDataProc = processBehaveLabels(behaveData, eventNameList)

for label_i = 1:length(eventNameList)
n_trials = size(behaveData, 3);
behaveDataProc.(eventNameList{label_i}) = zeros(n_trials, 1);
for k = 1:n_trials
    behaveDataProc.(eventNameList{label_i})(k) = any(any(behaveData(:,:,k) == label_i));    
end
end



histBehave = zeros(length(eventNameList)-4, size(behaveData, 1), size(behaveData, 3));

for t = 1:size(behaveData, 1)
    for label_i = 1:length(eventNameList)-4
        for T = 1:size(behaveData, 3)
            ind = find(strcmp(eventNameList, eventNameList{label_i}));
            histBehave(label_i, t, T) = histBehave(label_i, t, T) + sum(sum(behaveData(t, :, T) == ind));
        end
    end
end

behaveDataProc.trialBehave = zeros(size(behaveData, 1), size(behaveData, 3));
for t = 1:size(behaveData, 1)
    for T = 1:size(behaveData, 3)
        [maxVal, behaveDataProc.trialBehave(t, T)] = max(histBehave(:, t, T));
        if maxVal == 0
            behaveDataProc.trialBehave(t, T) = 0;
        end
    end
end


