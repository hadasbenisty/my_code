function plotNeuronsClustersWithMartkers(allMat_all)

M=[];
for k=1:length(allMat_all)
   M = [M; allMat_all{k}];
end
figure;imagesc(M);colormap gray;
L = 0;
for k=1:length(allMat_all)
    line([0 size(allMat_all{k}, 2)], [1 1]*L,'Color','b' );
    L = L + size(allMat_all{k}, 1);
end