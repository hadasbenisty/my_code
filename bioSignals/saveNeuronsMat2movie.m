function saveNeuronsMat2movie(behaveDataProc, mat2display, isbusy, timevec, roiLocs, avifilename)

frameReps=5;
N = length(mat2display);
indsmat=cell(N,1);
behavehist=cell(N,1);
for n=1:N
indsmat{n}=gray2ind(mat2gray(mat2display{n}));
behavehist{n} = hist(behaveDataProc{n}.trialBehave', 0:4);
    behavehist{n} = behavehist{n}([2 3 5],:);
end
figure;

for t_ind = 1:length(timevec)
    
    clustering = isbusy(:,t_ind);
    clusters = unique(clustering);
    resp=[];indsvec=[];
    alllocs = cell2mat(transpose(roiLocs{1}));
    clrmap=colormap(gray);
    for ci=1:length(clusters)
        ind = find(clustering==clusters(ci));
        for n=1:N
            resp{n}(ci) = mat2display{n}(ind(1),t_ind);
            indsvec{n}(ci) = indsmat{n}(ind(1), t_ind);
        end
    end
    for n=1:N
        subplot(2,N,n)
        imagesc(linspace(min(alllocs(:,1))-50, max(alllocs(:,1))+50, 100),linspace(min(alllocs(:,2))-50, max(alllocs(:,2))+50, 100), ones(100));colormap(jet);    hold all;
        hold all;
        for ci = 1:length(resp{n})
            celldata = transpose(roiLocs{1});
            matdata = cell2mat(celldata(clustering == clusters(ci)));
            plot(matdata(:,1),matdata(:,2),'LineStyle','none', 'Marker','o','MarkerEdgeColor',clrmap(indsvec{n}(ci)+1,:),...
                'MarkerFaceColor',clrmap(indsvec{n}(ci)+1,:),'MarkerSize',12)
            hold all;
            
            
        end
        hold off;
subplot(2,2*N, 2*N+2*n-1)
imagesc(timevec, 1:3,behavehist{n});set(gca, 'YTick',1:3); set(gca, 'YTickLabel',{'Lift','Grab','At Mouth'});
line([1 1]*4, [0 4],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
xlabel('Time[sec]');
    line([1 1]*timevec(t_ind), [0 4],'Color','g', 'lineWidth',2);
subplot(2,2*N, 2*N+2*n)
imagesc(timevec, 1:size(mat2display{n},1),mat2display{n});
line([1 1]*4, [0 size(mat2display{n},1)],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
xlabel('Time[sec]');ylabel('Neurons');
    line([1 1]*timevec(t_ind), [0 size(mat2display{n},1)],'Color','g', 'lineWidth',2);colormap gray

    end
    suptitle(sprintf('Location of Neurons Clusters & Behavioral Hostograms\n t = %2.2f secs', timevec(t_ind) ));
   Mtiling(t_ind) = getframe(gcf);
end
savemove2avi(avifilename, Mtiling, frameReps);