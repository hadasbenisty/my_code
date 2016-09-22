function [Trees, currSolutionTiling, neurons_order, time_order, dual_aff, meanDataLog] = runQuestAndTile(data, params, p1, runningOrder)

%% Ques.
[ Trees, dual_aff] = RunGenericDimsQuestionnaire( params, permute(data,(runningOrder) ) );


%% Tiling

[~, neurons_order] = sort(Trees{runningOrder==1}{2}.clustering);
time_order = getTimeOrderFromTree(Trees{runningOrder==2}{2}.clustering);
row_orderedtree = organizeTreeForTiling(Trees{runningOrder==1}, neurons_order);
col_orderedtree = organizeTreeForTiling(Trees{runningOrder==2}, time_order);
orderedData = data(neurons_order, :, :);%
orderedData = orderedData(:, time_order, :);
meanData = mean(orderedData, 3);
meanDataLog = log(threshold(meanData+2*abs(min(meanData(:))),0)+eps);



if params.verbose == 2
    figure;
    subplot(2,2,1);
    plotTreeWithColors(Trees{runningOrder==1}, 1:size(data,1));
    title('Neurons');
    subplot(2,2,2);
    plotTreeWithColors(Trees{runningOrder==2}, 1:size(data,2));
    title('Time');
    subplot(2,2,3);
    imagesc(meanDataLog);
    subplot(2,2,4);
end
[minCurrErr, currSolutionTiling] = reversedIterativeTilingGenericDimsleagal(meanDataLog, {row_orderedtree, col_orderedtree}, p1);
