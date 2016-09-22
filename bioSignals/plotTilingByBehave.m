function plotTilingByBehave(meanDataLog, X_N, eventNameList, behaveDataProc, roiLocs, neurons_order, time_order, solutionTiling, byttl, timevec, animalName, behave_events_times, allLabels,avifilename)



nr = length(roiLocs{1});

figure;
imagesc(timevec, 1:nr, meanDataLog);
line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
xlabel('Time[sec]');
ylabel('Neurons');colormap gray;
title(byttl);
for k=1:length(behave_events_times)
    line([1 1]*behave_events_times(k), [0 nr],'Color','g', 'lineWidth',2);
end
allbehavehist = hist(allLabels.trialBehave', 0:4);
allbehavehist = allbehavehist([2 3 5],:);

figure;
subplot(2,1,1);imagesc(timevec, 1:nr, meanDataLog);
line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
xlabel('Time[sec]');
ylabel('Neurons');colormap gray;
title('Samples')
subplot(2,1,2);
plot(timevec, allbehavehist);axis tight;
line([1 1]*4, [0 max(allbehavehist(:))],'Color','r', 'lineWidth',1);
title('Behavioral Data');
legend({'Lift','Grab','At Mouth'});
suptitle(byttl);

figure;
a(1) = subplot(2,1,1);
mat2display = getTiledData(meanDataLog, solutionTiling);

imagesc(timevec, 1:nr, mat2display);
title('Tiled Data'); colormap gray;
line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
for k=1:length(behave_events_times)
    line([1 1]*behave_events_times(k), [1 nr],'Color','g', 'lineWidth',3);
end
xlabel('Time[sec]');
ylabel('Neurons');

a(2) = subplot(2,1,2);

plot(timevec, allbehavehist);axis tight;
line([1 1]*4, [0 max(allbehavehist(:))],'Color','r', 'lineWidth',1);
for k=1:length(behave_events_times)
    line([1 1]*behave_events_times(k), [0 max(allbehavehist(:))],'Color','g', 'lineWidth',3);
end
title('Behavioral Data - All');

xlabel('Time[sec]');


suptitle([ animalName ' All Data, Tiled By ' byttl]);
linkaxes(a, 'x');legend({'Lift','Grab','At Mouth'});


for k = 1:length(X_N)
    orderedData_all{k} = X_N{k}(neurons_order, :, :);%#ok<*AGROW> %
    orderedData_all{k} = orderedData_all{k}(:, time_order, :);
    
    meanData_all{k}.alltrials = mean(orderedData_all{k}, 3);
    meanDataLog_all{k}.alltrials = log(threshold(meanData_all{k}.alltrials,0)+eps);
    for event_i = 1:length(eventNameList{k})
        if any(behaveDataProc{k}.(eventNameList{k}{event_i}))
            meanData_all{k}.(eventNameList{k}{event_i}) = mean(orderedData_all{k}(:,:,behaveDataProc{k}.(eventNameList{k}{event_i})==1), 3);
            meanDataLog_all{k}.(eventNameList{k}{event_i}) = log(threshold(meanData_all{k}.(eventNameList{k}{event_i})+2*abs(min(meanData_all{k}.(eventNameList{k}{event_i})(:))),0)+eps);
            
        end
    end
end

figure;
for k  = 1:length(meanDataLog_all)
    subplot(2,length(meanDataLog_all),k)
    mat2display = getTiledData(meanDataLog_all{k}.alltrials, solutionTiling);
    
    imagesc(timevec, 1:nr, mat2display);
    title(['Exp no. ' num2str(k)]); colormap gray;
    line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
    xlabel('Time[sec]');
    ylabel('Neurons');
    
    subplot(2,length(meanDataLog_all),length(meanDataLog_all)+k)
    %         imagesc(timevec, 1:size(behaveDataProc{k}.trialBehave,2), behaveDataProc{k}.trialBehave(toneTimeFramesExrtazoom,:)');
    %         title('Behave');
    imagesc(timevec, 1:nr, mat2display);
    title(['Exp no. ' num2str(k) ' Zoom']); colormap gray;
    line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        xlim([3.95 4.7]);xlabel('Time[sec]');
    ylabel('Neurons');
    
end
suptitle([ animalName ' All Data, Tiled By ' byttl]);
figure;
for k  = 1:length(meanDataLog_all)
    a(k) = subplot(2,length(meanDataLog_all),k);
    mat2display = getTiledData(meanDataLog_all{k}.alltrials, solutionTiling);
    
    imagesc(timevec, 1:nr, mat2display);
    title(['Exp no. ' num2str(k)]); colormap gray;
    line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
    xlabel('Time[sec]');
    ylabel('Neurons');
    
    a(length(meanDataLog_all)+k) = subplot(2,length(meanDataLog_all),length(meanDataLog_all)+k);
    
    plot(timevec, allbehavehist);axis tight;
    line([1 1]*4, [0 max(allbehavehist(:))],'Color','r', 'lineWidth',1);
    title('Behavioral Data - All');
    
    xlabel('Time[sec]');
    
end
suptitle([ animalName ' All Data, Tiled By ' byttl]);
linkaxes(a, 'x');legend({'Lift','Grab','At Mouth'});
if 0
    for event_i = 1:length(eventNameList{k})
        f=figure;
        for k  = 1:4
            if any(behaveDataProc{k}.(eventNameList{k}{event_i}))
                subplot(2,4,k);
                
                mat2display = getTiledData(meanDataLog_all{k}.(eventNameList{k}{event_i}), solutionTiling);
                
                imagesc(timevec, 1:nr, mat2display);
                title(['Exp no. ' num2str(k)]); colormap gray;
                line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
                
                subplot(2,4,4+k)
                %                 imagesc(timevec, 1:sum(behaveDataProc{k}.(eventNameList{k}{event_i})), behaveDataProc{k}.trialBehave(toneTimeFramesExrtazoom,behaveDataProc{k}.(eventNameList{k}{event_i})==1)');
                %                 title('Behave');
                imagesc(timevec, 1:nr, (mat2display));
                title(['Exp no. ' num2str(k) ' Zoom']); colormap gray;
                line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);
                
                xlim([3.95 4.7]);
                issupttl=true;
            else
                close(f);
                issupttl=false;
                break;
            end
        end
        if issupttl
            suptitle([ animalName ' ' eventNameList{k}{event_i} ' Data, Tiled By ' byttl]);
        end
    end
end
figure;
for k  = 1:length(meanDataLog_all)
    behavehist = hist(behaveDataProc{k}.trialBehave', 0:4);
    behavehist = behavehist([2 3 5],:);
    
    a(k) = subplot(2,length(meanDataLog_all),k);
    mat2display = getTiledData(meanDataLog_all{k}.alltrials, solutionTiling);
    
    imagesc(timevec, 1:nr, mat2display);
    title(['Exp no. ' num2str(k)]); colormap gray;
    line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
    xlabel('Time[sec]');
    ylabel('Neurons');
    a(k+length(meanDataLog_all)) = subplot(2,length(meanDataLog_all),k+length(meanDataLog_all));
    plot(timevec, behavehist);axis tight;
    line([1 1]*4, [0 max(behavehist(:))],'Color','r', 'lineWidth',1);
end
legend({'Lift','Grab','At Mouth'});
linkaxes(a, 'x');

figure;
for k  = 1:4
    
    subplot(2,2,k)
    mat2display = getTiledData(meanDataLog_all{k}.alltrials, solutionTiling);
    
    imagesc(timevec, 1:nr, mat2display);
    title(['Exp no. ' num2str(k)]); colormap gray;
    line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
    xlabel('Time[sec]');
    ylabel('Neurons');
    for l=1:length(behave_events_times)
        line([1 1]*behave_events_times(l), [0 nr],'Color','g', 'lineWidth',2);
    end
end



suptitle([ animalName ' All Data, Tiled By ' byttl]);
%% make a movie of the neurons
if ~isempty(avifilename)
    for k=1:4
        tiledData{k} = getTiledData(meanDataLog_all{k}.alltrials, solutionTiling);
        rawdata{k} =meanDataLog_all{k}.alltrials;
    end
    for k=1:4
        saveNeuronsMat2movie(behaveDataProc(k), tiledData(k), solutionTiling.isbusy, timevec, roiLocs, [avifilename, 'tiling_exp_no' num2str(k) '.avi']);
    end
    saveNeuronsMat2movie(behaveDataProc, tiledData, solutionTiling.isbusy, timevec, roiLocs, [avifilename, 'tiling.avi']);
    
    figure
    for k=1:4
        subplot(2,2,k)
        mat2display = getTiledData(meanDataLog_all{k}.(eventNameList{k}{event_i}), solutionTiling);
        
        imagesc(timevec, 1:nr, log(threshold(mat2display,.6)));xlabel('Time[sec]');
        ylabel('Neurons');
        title(['Exp no. ' num2str(k) ' Zoom']); colormap gray;
        line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        xlim([3.95 4.2]);
    end
    figure;
    for k=1:4
        subplot(2,2,k)
        mat2display = getTiledData(meanDataLog_all{k}.(eventNameList{k}{event_i}), solutionTiling);
        
        imagesc(timevec, 1:nr, log(threshold(mat2display,.6)));
        title(['Exp no. ' num2str(k) ' Zoom']); colormap gray;
        line([1 1]*4, [1 nr],'Color','r', 'lineWidth',2);        xlim([4.8 6]);
    end
end