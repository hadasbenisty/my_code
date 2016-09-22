function plotHistByLevel(tree, level, labels)


clusters = unique(tree{level}.clustering);
R = ceil(sqrt(length(clusters)));
for ci=1:length(clusters)
subplot(R,R,ci);

c = categorical(labels(tree{level}.clustering == clusters(ci)));
histogram(c);title(['Cluster no. ' num2str(ci)]);
end
suptitle(['Tree Level ' num2str(level)]);
end
