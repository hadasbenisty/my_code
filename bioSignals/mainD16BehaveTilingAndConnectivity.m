%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script for analyzing the samples of D16 - quest, tiling, behavioral 
% 
%
%
% Written by Hadas Benisty, 9/22/2016
% hadas.benisty@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function mainD30SucFailNewSamples(dims, epsilonInd, overwrite, runOnlyOnToneFrames)
clear ;
% ca;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../../Matlab_Utils/'));
addpath(genpath('../tiling'));
addpath(genpath('../../grangerCausality/'));
addpath(genpath('../../SvdKmeansCluster'));
addpath(genpath('../../GraphConnectivity/'));
%% D30

animalName = 'D16';
datapth = ['../../../datasets/biomed/' animalName ];
runOnlyOnToneFrames = 2;epsilonInd=9;overwrite=0;dims=2;

isPermNeurons = true;

nt = 360;
files = {'6_24_14_1-15','7_5_14_1-25','7_9_14_1-18','7_15_15_1-46','7_17_14_1-25','7_22_14_1-16' '8_5_14_1-51'  };%
toneVecAlltime = [ones(1,120) 100*ones(1, 240)];



%% Set Params
% selectedTimeFrams = 115:140;
dims4quest=3;
params = SetGenericDimsQuestParams(dims4quest, true);
for ind = 1:dims4quest
    params.init_aff{ind}.metric = 'cosine_similarity';
    params.tree{ind}.runOnEmbdding = ~false;
    params.tree{ind}.splitsNum=2;
    params.tree{ind}.treeDepth = 8;%4;

end
runningOrder = [  3 2 1];

params.tree{runningOrder==2}.eigs_num=50;
params.tree{runningOrder==3}.eigs_num=50;

params.tree{runningOrder==2}.splitsNum = [8 2 2 2 2 2 2 2 2];%9
params.tree{runningOrder==3}.runOnEmbdding = true;
params.n_iters = 2;
params.tree{runningOrder==2}.treeDepth = 5;
params.data.over_rows = false;
params.data.to_normalize = true;%false
params.data.normalization_type = 'by_std';
params.n_iters=2;

t = linspace(0, 12, nt);
% toneTimeFramesExrtazoom = 115:160;%115:140;
toneTimeFramesExrtazoom =100:185;
%% Load Data
X = cell(length(files), 1); NeuronsLabels = cell(length(files), 1);

% [X_all, expLabel_all, NeuronsLabels_all, behaveDataProcessed, bahavetitles] = loadNeuronsData(datapth, files, nt);
% XalltimeAll = X_all;

for n = 1:length(files)
    s = load(fullfile(datapth, files{n}, 'sampsAndBehave'));
    Xallneurons{n} = s.dffDataArray;
    behaveData{n} = s.eventDataArray(toneTimeFramesExrtazoom, :,:);
    NeuronsLabels{n} = s.roiNames;
    eventNameList{n} = s.eventNameList;
    roiLocs{n} = s.roiLoc;
end
combLabels = NeuronsLabels{1};
for n=2:length(files)
    combLabels = intersect(combLabels, NeuronsLabels{n});
end
% rand_permuation = randperm(length(combLabels));
rand_permuation = 1:length(combLabels);

for n=1:length(files)
    combinds = getCombInds(NeuronsLabels{n}, combLabels);    
    Xalltime{n} = Xallneurons{n}(combinds, :, :);
    if isPermNeurons
       Xalltime_perm{n} = Xalltime{n}(rand_permuation, :, :);
    else
        Xalltime_perm{n} = Xalltime{n};
    end
    NeuronsLabels_all{n} = NeuronsLabels{n}(combinds);
    roiLocs_all{n} = roiLocs{n}(combinds);
end

toneVec =  toneVecAlltime(toneTimeFramesExrtazoom);
for n = 1:length(files)
    X{n} = Xalltime_perm{n}(:, toneTimeFramesExrtazoom, :);
end

X_N_all=[];
X_N = cell(length(files), 1);
for n = 1:length(X)
    X_N{n} = normalizeDataWithOrdering(X{n}, runningOrder, params.data);
    
end
for n = 1:length(X)
X_N_all = cat(3,X_N_all,X_N{n});
end
X_N_all = normalizeDataWithOrdering(X_N_all, runningOrder, params.data);
en = 0;
for n = 1:length(X)
    st = en +1 ;
    en = en + size(X_N{n},3);
    X_NN{n} = X_N_all(:,:,st:en);
end
params.data.to_normalize = false;



%% process behaveData{n} into labels
for n = 1:length(X)
    for label_i = 1:length(eventNameList{n})
        allLabels.(eventNameList{n}{label_i})=[];
    end
end
allLabels.trialBehave=[];
for n = 1:length(behaveData)
    behaveDataProc{n} = processBehaveLabels(behaveData{n}, eventNameList{n});
end
for n = 1:4
    for label_i = 1:length(eventNameList{n})
        allLabels.(eventNameList{n}{label_i}) = [allLabels.(eventNameList{n}{label_i}); behaveDataProc{n}.(eventNameList{n}{label_i})];
    end
    allLabels.trialBehave = [allLabels.trialBehave behaveDataProc{n}.trialBehave];
end
nr = size(X{1},1);
behave_events_times_all = getBehaveTImeAxis(t(toneTimeFramesExrtazoom), allLabels.trialBehave, 2);
% figure;
% imagesc(t, 1:size(allLabels.trialBehave,2),allLabels.trialBehave');colormap gray;
% title('Behavioral Data');
% line([1 1]*4, [1 size(allLabels.trialBehave,2)],'Color','r', 'lineWidth',2);        %xlim([3.95 4.7]);
% xlabel('Time[sec]');
% ylabel('Trials');

% figure;
% imagesc(t(toneTimeFramesExrtazoom), 1:size(allLabels.trialBehave,2),allLabels.trialBehave');colormap gray;
% title('Behavioral Data');
% line([1 1]*4, [1 size(allLabels.trialBehave,2)],'Color','r', 'lineWidth',2);       colormap gray;
% xlabel('Time[sec]');
% ylabel('Trials');
% for k=1:length(behave_events_times_all)
% line([1 1]*behave_events_times_all(k), [0 size(allLabels.trialBehave,2)],'Color','g', 'lineWidth',2);        xlim([3.95 6]);
% end

p1.representation_error_fun = @squaredErr;
p1.verbose = 0;params.verbose = 0;


eps = [0.05 0.15 0.1 0.08 0.05 0.1 0.07];

%% Run Qu. + Tiling
% quest and tile using all data
p1.eps = 0.1;
runningOrder = [ 1 2 3];
[ Trees, dual_aff] = RunGenericDimsQuestionnaire( params, permute(X_N_all,(runningOrder) ) );
[vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==3}, 0.0), 3);
embedding{runningOrder==3} = vecs*vals;

PlotEmbedding(embedding{runningOrder==3}, [ones(14,1); 2*ones(25,1); 3*ones(18,1); 4*ones(46, 1); 5*ones(25,1); 6*ones(16,1); 7*ones(50,1)],  ['Dim No. ' num2str(2) ' - Final '] );

[Trees_all, currSolutionTiling_all, neurons_order_all, time_order_all, ...
    dual_aff_all, meanDataLog_all] = runQuestAndTile(X_N_all, params, p1, runningOrder);
PlotEmbedding(embedding{runningOrder==3}, [ones(14,1); ones(25,1);ones(18,1); ones(46, 1); ones(25,1); 2*ones(16,1); 2*ones(50,1)],  ['Dim No. ' num2str(2) ' - Final '] );

allbehavehist = hist(allLabels.trialBehave', 0:4);
allbehavehist = allbehavehist([2 3 5],:);


for n = length(X_NN):-1:1
    p1.eps=eps(n);
    [Trees{n}, currSolutionTiling{n}, neurons_order{n}, time_order{n}, dual_aff{n}, meanDataLog{n}] = ...
        runQuestAndTile(X_N{n}, params, p1, runningOrder);
behave_events_times{n} = getBehaveTImeAxis(t(toneTimeFramesExrtazoom), behaveDataProc{n}.trialBehave, 2);

    plotTilingByBehaveD16(meanDataLog{n}, X_NN, eventNameList, behaveDataProc, roiLocs_all, ...
    neurons_order{n}, time_order{n}, currSolutionTiling{n}, num2str(n), t(toneTimeFramesExrtazoom), animalName, behave_events_times{n}, behaveDataProc{n},'');
end
figure;
for k = 1:length(X_NN)
    suc(k) = sum(behaveDataProc{k}.success)/length(behaveDataProc{k}.success);
    fail(k) = sum(behaveDataProc{k}.failure)/length(behaveDataProc{k}.failure);

    subplot(4,2,k);
    
    behavehist = hist(behaveDataProc{k}.trialBehave', 0:14);
    disp(unique(behaveDataProc{k}.trialBehave)');
behavehist = behavehist(2:end,:);
plot(behavehist');
title(num2str(k));
end
legend(eventNameList{n});%, 'Location','BestOutside'); % (do not touch));
subplot(4,2,8);
plot(suc);
xlabel('Exp');
ylabel('Suc. Rate');

meanTrials = false;
params.connectivity.winLen = 1;
params.connectivity.overlap = .5;
params.connectivity.epsParamD = 1;
params.connectivity.epsParamVis = 1;
params.connectivity.distTypeAlltimes = 3;
params.connectivity.distTypeTimeSample = 4;
params.connectivity.meanTrials = false;
params.scattering.winlen = 0;
permutedData = mean(X_N_all,3);
meanpermutedDataLog = log(threshold(permutedData+2*abs(min(permutedData(:))),0)+eps);
permutedDataLog = permute(log(threshold(X_N_all+2*abs(min(X_N_all(:))),0)+eps),[2 1 3]);

S_table_permute = scatteringWrapper(meanpermutedDataLog, params.scattering.winlen);
input_scattering_permuted = permute(S_table_permute, [  1 3 2]);
if meanTrials
    [ DallTimes_permuted, D_neurons_permuted, Vis_neurons_permuted, VisNorm_neurons_permuted, Vis_Dalltimes_neurons_permuted ] = ...
    GraphConnectivity( meanpermutedDataLog.',  params.connectivity);
else
    [ DallTimes_permuted, D_neurons_permuted, Vis_neurons_permuted, VisNorm_neurons_permuted, Vis_Dalltimes_neurons_permuted ] = ...
    GraphConnectivity( permutedDataLog,  params.connectivity);

% for trial_i=1:size(permutedDataLog, 3)
%      [ DallTimes_permuted_trials{trial_i}, D_neurons_permuted_trials{trial_i}, Vis_neurons_permuted_trials{trial_i}, VisNorm_neurons_permuted_trials{trial_i}...
%          , Vis_Dalltimes_neurons_permuted_trials{trial_i} ] = GraphConnectivity( permutedDataLog(:,:,trial_i),  params.connectivity);
% 
% end
end
    % [ DallTimes_permuted, D_neurons_permuted, Vis_neurons_permuted, VisNorm_neurons_permuted, Vis_Dalltimes_neurons_permuted ] = ...
%     GraphConnectivity( input_scattering_permuted,  params.connectivity);

DallTimes_neurons = DallTimes_permuted(neurons_order_all, neurons_order_all);
D_neurons = D_neurons_permuted(:,:,neurons_order_all );
D_neurons = D_neurons(:,neurons_order_all,: );
Vis_neurons = Vis_neurons_permuted(:,:, neurons_order_all);
Vis_neurons = Vis_neurons(:,neurons_order_all, :);
Vis_Dalltimes_neurons = Vis_Dalltimes_neurons_permuted(neurons_order_all, neurons_order_all);

% for trial_i=1:size(permutedDataLog, 3)
%    DallTimes_neurons_trials(:,:,trial_i) = DallTimes_permuted_trials{trial_i}(neurons_order_all, neurons_order_all);
% D_neurons_trials(:,:,:,trial_i) = D_neurons_permuted_trials{trial_i}(:,:,neurons_order_all );
% D_neurons_trials(:,:,:,trial_i) = D_neurons_trials(:,neurons_order_all,:,trial_i);
% Vis_neurons_trials(:,:,:,trial_i) = Vis_neurons_permuted_trials{trial_i}(:,:, neurons_order_all);
% Vis_neurons_trials(:,:,:,trial_i) = Vis_neurons_trials(:,neurons_order_all,:,trial_i);
% Vis_Dalltimes_neurons_trials(:,:,trial_i) = Vis_Dalltimes_neurons_permuted_trials{trial_i}(neurons_order_all, neurons_order_all);
% end

figure;
subplot(2,1,1);
imagesc(Vis_Dalltimes_neurons_permuted); colormap('gray');
title('Unsorted'); xlabel('Neurons'); ylabel('Neurons');
subplot(2,1,2);
imagesc(Vis_Dalltimes_neurons); colormap('gray');
title('Sorted By Trees'); xlabel('Neurons'); ylabel('Neurons')
 mysave(gcf,'neurons_conn')
params.connectivity.distTypeAlltimes=3;
params.connectivity.distTypeTimeSample = 4;
params.connectivity.winLen = 1;
params.connectivity.overlap = .5;
if meanTrials
    [ DallTimes_time, D_time, Vis_time, VisNorm_time, Vis_Dalltimes_time ] = ...
    GraphConnectivity( meanpermutedDataLog,  params.connectivity);

else
[ DallTimes_time, D_time, Vis_time, VisNorm_time, Vis_Dalltimes_time ] = ...
    GraphConnectivity( permute(permutedDataLog, [2 1 3]),  params.connectivity);
end
% [ DallTimes_time, D_time, Vis_time, VisNorm_time, Vis_Dalltimes_time ] = ...
%     GraphConnectivity( permute(input_scattering_permuted, [2 1 3]),  params.connectivity);

Vis_time = Vis_time(neurons_order_all,:,:);

[~, neurons_order] = sort(Trees_all{runningOrder==1}{2}.clustering);
time_order = getTimeOrderFromTree(Trees_all{runningOrder==2}{2}.clustering);
row_orderedtree = organizeTreeForTiling(Trees_all{runningOrder==1}, neurons_order);
col_orderedtree = organizeTreeForTiling(Trees_all{runningOrder==2}, time_order);
[Nr, T] = size(meanpermutedDataLog);
figure;subplot(2,2,4);
imagesc(t(toneTimeFramesExrtazoom), t(toneTimeFramesExrtazoom), Vis_Dalltimes_time); xlabel('time');
hold all;

trans_t = find(diff(col_orderedtree{2}.clustering));
for n_i = 1:length(trans_t)
    line(t(toneTimeFramesExrtazoom(trans_t(n_i)))*[1 1], [0 T],'Color','b' );
    line([0 T], t(toneTimeFramesExrtazoom(trans_t(n_i)))*[1 1], 'Color','b' );
end
title('Conncetivity  - Time'); xlabel('time'); ylabel('time');
 subplot(2,2,2); plot(t(toneTimeFramesExrtazoom), allbehavehist');axis tight;title('Behaviroal Histogram');legend('Lift','Grab','At Mouth');
subplot(2,2,3);
imagesc(Vis_Dalltimes_neurons); xlabel('time');
hold all;
trans_nr = find(diff(row_orderedtree{2}.clustering));
for n_i = 1:length(trans_nr)
    line(trans_nr(n_i)*[1 1], [0 Nr],'Color','b' );
    line([0 Nr], trans_nr(n_i)*[1 1], 'Color','b' );
end
title('Conncetivity  - Neurons'); xlabel('neurons'); ylabel('neurons');
subplot(2,2,1); treeplot(nodes(row_orderedtree), '.');
title('Neurons Tree');
if meanTrials
 mysave(gcf,'trials_meanBeforeEverything');
else
     mysave(gcf,'trials_kernel');
end

% figure;   
%  for n_i=1:size(Vis_neurons, 1)     
%     subplot(2,1,1);
%      imagesc(squeeze(real(Vis_neurons(n_i,:,:))));
%     subplot(2,1,2);
%    plot(t(toneTimeFramesExrtazoom), allbehavehist');axis tight
%    hold all;
%    line(t(toneTimeFramesExrtazoom(n_i))*[1 1], [0 max(allbehavehist(:))], 'Color','g');
%      hold off;
%     d(t_i)=getframe(gcf);
% %     M(t_i) = im2frame(d.cdata,d.colormap);
%  end
%  savemove2avi('hadas', d, 1);
%  
 tiledData = getTiledData(meanDataLog_all, currSolutionTiling_all);

  figure; 
  Nr = size(currSolutionTiling_all.isbusy,1);
  time_transitions = [1 find(diff(currSolutionTiling_all.isbusy(1,:))) size(Vis_neurons, 1)] ;
 for t_i=1:length(time_transitions)-1
     trans_nr = find(diff(currSolutionTiling_all.isbusy(:,time_transitions(t_i)+1))~=0);
     midtimetrans = mean(t(toneTimeFramesExrtazoom([time_transitions(t_i),time_transitions(t_i+1)])));
     for t_reps = time_transitions(t_i):time_transitions(t_i+1)
     
     % data 
     subplot(2,2,1);
      imagesc(t(toneTimeFramesExrtazoom),  1:Nr, meanDataLog_all);
 hold all;
 line([1 1]*4, [0 Nr],'Color','r', 'lineWidth',1);xlabel('Time[secs]');title('Data');ylabel('Neurons');
line(t(toneTimeFramesExrtazoom(t_reps))*[1 1], [0 Nr], 'Color','g');
% tiled data
  subplot(2,2,2);
      imagesc(t(toneTimeFramesExrtazoom),  1:Nr, tiledData);
 hold all;
 line([1 1]*4, [0 Nr],'Color','r', 'lineWidth',1);xlabel('Time[secs]');ylabel('Neurons');title('Tiled Data');
line(t(toneTimeFramesExrtazoom(t_reps))*[1 1], [0 Nr], 'Color','g');

   
     % behaviroal
   subplot(2,2,3);
   
     
   plot(t(toneTimeFramesExrtazoom), allbehavehist');axis tight
   hold all;
   line(t(toneTimeFramesExrtazoom(t_reps))*[1 1], [0 max(allbehavehist(:))], 'Color','g');
     hold off;xlabel('Time[secs]');
     % connectivity
     subplot(2,2,4);
       imagesc(squeeze(mean(real(Vis_neurons(time_transitions(t_i):time_transitions(t_i+1),:,:)),1)));colormap gray
       hold on;
     for trans_i = 1:length(trans_nr)
         line(trans_nr(trans_i)*[1 1], [0 Nr], 'Color','b', 'lineWidth',1);
         line([0 Nr], trans_nr(trans_i)*[1 1],  'Color','b', 'lineWidth',1);
     end
     hold off;
        xlabel('Neurons');
     ylabel('Neurons');
     title('Connectivity');
    d(t_i)=getframe(gcf);
     end
%     M(t_i) = im2frame(d.cdata,d.colormap);
 end
 savemove2avi('tiling', d, 1);
% 
%  
%   Nr = size(currSolutionTiling_all.isbusy,1);
%  for t_i=[10 19:38 29 34 62]
%      figure
%      trans_nr = find(diff(currSolutionTiling_all.isbusy(:,t_i))~=0);
%     subplot(2,2,1);
%       imagesc(t(toneTimeFramesExrtazoom),  1:Nr, meanDataLog_all);
%  hold all;
%  line([1 1]*4, [0 Nr],'Color','r', 'lineWidth',1);xlabel('Time[secs]');title('Data');ylabel('Neurons');
% line(t(toneTimeFramesExrtazoom(t_i))*[1 1], [0 Nr], 'Color','g');
% 
%      subplot(2,2,2);
%       imagesc(t(toneTimeFramesExrtazoom),  1:Nr, tiledData);
%  hold all;
%  line([1 1]*4, [0 Nr],'Color','r', 'lineWidth',1);xlabel('Time[secs]');ylabel('Neurons');title('Tiled Data');
% line(t(toneTimeFramesExrtazoom(t_i))*[1 1], [0 Nr], 'Color','g');
% 
%      
%     subplot(2,2,3);
%    
%      
%    plot(t(toneTimeFramesExrtazoom), allbehavehist');axis tight
%    hold all;
%    line(t(toneTimeFramesExrtazoom(t_i))*[1 1], [0 max(allbehavehist(:))], 'Color','g');
%      hold off;xlabel('Time[secs]');
%      subplot(2,2,4);
%       imagesc(squeeze(real(Vis_neurons(t_i,:,:))));colormap gray
%      for trans_i = 1:length(trans_nr)
%          line(trans_nr(trans_i)*[1 1], [0 Nr], 'Color','b', 'lineWidth',1);
%          line([0 Nr], trans_nr(trans_i)*[1 1],  'Color','b', 'lineWidth',1);
%      end
%      xlabel('Neurons');
%      ylabel('Neurons');
%      title('Connectivity');
%    suptitle(['t = ' num2str(t(toneTimeFramesExrtazoom(t_i))) '[secs]']);
%  end
%  
%  tilinginds = [1; find(diff(currSolutionTiling_all.isbusy(:,24))); Nr];
%  for nr=1:length(tilinginds)-1
%  rowNum = tilinginds(nr):tilinginds(nr+1)-1;
% figure
% imagesc(t(toneTimeFramesExrtazoom),1:97,squeeze(mean(real(D_neurons(:,rowNum,:)),2))); colormap('gray')
% xlabel('time'); ylabel('Neurons'); title(['Conncetivity over time between all neurons in tile #',num2str(nr)])
%  end
%  
%  
%  
%  %%
%  
%   figure; 
%   T = size(currSolutionTiling_all.isbusy,2);
%   neurons_transitions = [find(diff(currSolutionTiling_all.isbusy(:,1))); size(Vis_time, 1)] ;
%  for nr=1:length(neurons_transitions)-1
%      trans_time = find(diff(currSolutionTiling_all.isbusy(neurons_transitions(nr)+1,:))~=0);
%      midnrtrans = mean([neurons_transitions(nr),neurons_transitions(nr+1)]);
%      for nr_reps = neurons_transitions(nr):neurons_transitions(nr+1)
%      
%      % data 
%      subplot(2,2,1);
%       imagesc(t(toneTimeFramesExrtazoom),  1:Nr, meanDataLog_all);
%  hold all;
%  line([1 1]*4, [0 Nr],'Color','r', 'lineWidth',1);xlabel('Time[secs]');title('Data');ylabel('Neurons');
% line([0 T], nr_reps*[1 1], 'Color','g');
% % tiled data
%   subplot(2,2,2);
%       imagesc(t(toneTimeFramesExrtazoom),  1:Nr, tiledData);
%  hold all;
%  line([1 1]*4, [0 Nr],'Color','r', 'lineWidth',1);xlabel('Time[secs]');ylabel('Neurons');title('Tiled Data');
% line([0 T], nr_reps*[1 1], 'Color','g');
% 
%    
%      % behaviroal
%    subplot(2,2,3);
%    
%      
%    plot(t(toneTimeFramesExrtazoom), allbehavehist');axis tight
%  xlabel('Time[secs]');
%      % connectivity
%      subplot(2,2,4);
%        imagesc(squeeze(mean(real(Vis_time(neurons_transitions(nr):neurons_transitions(nr+1),:,:)),1)));colormap gray
%        hold on;
%      for trans_i = 1:length(trans_time)
%          line([0 T],trans_time(trans_i)*[1 1], 'Color','b', 'lineWidth',.5);
%          line(trans_time(trans_i)*[1 1], [0 T], 'Color','b', 'lineWidth',.5);
%      end
%      hold off;
%         xlabel('Neurons');
%      ylabel('Neurons');
%      title('Connectivity');
%     dt(nr)=getframe(gcf);
%      end
% %     M(t_i) = im2frame(d.cdata,d.colormap);
%  end
%  savemove2avi('tilingt', dt, 1);
