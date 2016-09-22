function tilingSol = setTilingByTreeLevel(Trees)
for n = 1:length(Trees)
    clusters{n} = unique(Trees{n}{1}.clustering);
end
tilingSol.isbusy = zeros(4);
tilingInd = 1;
cmdStr = '';
for n = 1:length(clusters)
    cmdStr = [cmdStr ' for k' num2str(n) ' = 1:length(clusters{' num2str(n) '}),'];
end
for n = 1:length(clusters)
    cmdStr = [cmdStr ' inds' num2str(n) '= Trees{' num2str(n) '}{1}.clustering == clusters{' num2str(n) '}(k' num2str(n) ');'];
end
cmdStr = [cmdStr 'tilingSol.isbusy('];
for n = 1:length(clusters)
cmdStr = [cmdStr 'inds' num2str(n) ','];
end
cmdStr = [cmdStr(1:end-1) ') = tilingInd; tilingInd = tilingInd + 1;'];
for n = 1:length(clusters)
cmdStr = [cmdStr 'end;'];
end
eval(cmdStr);

