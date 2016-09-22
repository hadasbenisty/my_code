function histByClustering(clustering, labels)


clusters = unique(clustering);
R = ceil(sqrt(length(clusters)));
for ci=1:length(clusters)
subplot(R,R,ci);

c = categorical(labels(clustering == clusters(ci)));
histogram(c);title(['Cluster no. ' num2str(ci)]);
end
end
