close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../Questionnaire'));
addpath(genpath('../utils'));
overwrite = false;
files = {'8_15_13_1-35' '8_12_14_1-40'};% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'
figspath1 = fullfile('D:\workWithBoss\summaries\D30',files{1});
for n=2:length(files)
    figspath1 = [figspath1 files{n}];
end
wrkspname = fullfile(figspath1, 'wrkspace.mat');
if exist(wrkspname, 'file') && ~overwrite
    load(wrkspname);
else
    rng(73631);
    %% Init params
    dorandperm_col = false;
    dorandperm_row = false;
    dorandperm_trials = false;
    eigsnum_col = 12;
    eigsnum_row = 12;
    eigsnum_trial= 12;
    row_alpha = .2;
    row_beta = 0;
    col_alpha = .2;
    col_beta = 0;
    trial_alpha = .2;
    trial_beta = 0;
    params  = SetQuest3DParams(eigsnum_col, eigsnum_row, eigsnum_trial, row_alpha, row_beta, col_alpha, col_beta, trial_alpha, trial_beta );
    % params = getParamsForBioData(eigsnum_col,eigsnum_row, eigsnum_trial);
    %% Load Data
    
    datapth = '..\..\..\datasets\biomed\D30';
    
    nt = 360;
    [X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);
    
    
    [nr, nt, nT] = size(X);
    % perm data
    if dorandperm_col
        col_perm = randperm(nt);
    else
        col_perm = 1:nt;
    end
    if dorandperm_row
        row_perm = randperm(nr);
    else
        row_perm = 1:nr;
    end
    if dorandperm_trials
        trial_perm = randperm(nT);
    else
        trial_perm = 1:nT;
    end
    data = X(row_perm, :, :);
    data = data(:, col_perm, :);
    data = data(:, :, trial_perm);
    
    %% Run Qu. 2D
    % data2D = permute(mean(data(:,:,1:4),3), [2 1 3]);
    % [ row_tree, col_tree ] = RunQuestionnaire( params, data2D );
    % row_aff = CalcEmdAff(data2D.', col_tree, params.row_emd);
    % col_aff = CalcEmdAff(data2D, row_tree, params.col_emd);
    %
    % [ err_rate ] = OrganizeData( data2D, data2D, row_aff, col_aff, col_perm, row_perm, 12, 12  );
    
    %% Run Qu. 3D
    [row_tree, col_tree, trials_tree, row_aff, col_aff, trials_aff] = RunQuestionnaire3D(params, data);
    % [col_tree,row_tree, trials_tree, col_aff, row_aff, trials_aff] = RunQuestionnaire3D(params, permute(data, [2 1 3]));
    % eigsnum_row = 3;
    % [ err_rate, organized_data, row_order, col_order, trial_order ] = OrganizeData3D( data, data, row_aff, col_aff, trials_aff, row_perm, col_perm, trial_perm, eigsnum_col, eigsnum_row, eigsnum_trial );
    
    %  coefs  = FindTreeAverages3D(data, col_tree, row_tree );
    
    %% evaluating cost at every node
    % X = data(:,:,1);
    
    %
    %
    % phi = 0.5;
    % W = exp(-phi * squareform(pdist(X.')));
    % W=W/sum(W(:))/sqrt(nt);
    % W_tild = exp(-phi * squareform(pdist(X)));
    % W_tild=W_tild/sum(W_tild(:))/sqrt(nr);
    %
    % gam = 100;
    % U = findLowestLevelAveragesByTrees(col_tree, row_tree, X);
    % cost = 0.5*sum(sum((X-U).^2)) + gam * (omega_fun(W, U) + omega_fun(W_tild, U.'));
    % % prune
    % col_tree1{1} = col_tree{1};
    % col_tree1{1}.clustering(col_tree{2}.clustering==1) = 1;
    % col_tree1{1}.folder_count = length(unique(col_tree1{1}.clustering));
    % U = findLowestLevelAveragesByTrees(col_tree1, row_tree, X);
    % cost = 0.5*sum(sum((X-U).^2)) + gam * (omega_fun(W, U) + omega_fun(W_tild, U.'));
    
    %% Visualization
    % organized_data1 = data(row_order,:, :);
    % organized_data1 = organized_data1(:,:, :);
    %
    % figure;
    % subplot(2,2,1);imagesc(shiftdim(data(:,:,1)))
    % subplot(2,2,2);imagesc(shiftdim(organized_data1(:,:,1)))
    % subplot(2,2,3);imagesc(shiftdim(data(:,:,2)))
    % subplot(2,2,4);imagesc(shiftdim(organized_data1(:,:,2)))
    % figure;subplot(3,1,1);
    % imagesc(shiftdim(mean(data(:,:,1:10), 3)))
    % subplot(3,1,2);
    % imagesc(shiftdim(mean(organized_data(:,:,1:10), 3)))
    % subplot(3,1,3);
    % imagesc(shiftdim(mean(organized_data1(:,:,1:10), 3)))
    % [row_vecs, row_vals] = CalcEigs(row_aff, 4);
    % figure;
    %
    % embedding = row_vecs*row_vals;
    %
    %
    % subplot(2,1,1);
    % PlotEmbedding( embedding, row_perm, ['Row' ] );
    % subplot(2,1,2);
    % plotEmbeddingWithColors(embedding, row_tree{4}.clustering, ['Row' ]);
    
    
    
    mkNewFolder(figspath1);
    savefigs = true;
    save(wrkspname);
end
% plotSlices(data, 9, 'Neurons','Time', 'Trials')
row_thresh = 0.0;
col_thresh = 0.0;
trials_thresh = 0;% 0.4
eigsnum_row = 3;
eigsnum_col = 3;
eigsnum_trials = 3;
tonetimeS = 100;
tonetimeE = 120;
toneLabel = [ones(tonetimeS, 1); 2*ones(tonetimeE-tonetimeS, 1); 3*ones(size(col_aff, 1)-(tonetimeE), 1)];
savefigs=true;
if length(unique(expLabel))==1
    expLabel=[];
end
% plotTreesAndEmbedding3D(figspath1, savefigs, 'Neurons', 'Time', 'Trials', ...
%     eigsnum_row, eigsnum_col, eigsnum_trials, ...
%     row_aff, col_aff, trials_aff, ...
%     row_thresh, col_thresh, trials_thresh, ...
%     row_tree, col_tree, trials_tree, ...
%     row_perm, col_perm, trial_perm,  [], toneLabel, expLabel);
% clc;
% for clusterLevel = 2:length(row_tree)-1
%     disp(['Tree Level ' num2str(clusterLevel)]);
%     levelclustering = row_tree{clusterLevel}.clustering;
%     clusters = unique(levelclustering);
%     clusteredData = cell(length(clusters), 1);
%     for ci = 1:length(clusters)
%         clusteredData{ci} = data(levelclustering==clusters(ci), :, :);
%         clustersinds{ci} = find(levelclustering==clusters(ci));
%         N = sum(levelclustering==clusters(ci));
%         R = ceil(sqrt(N));
%         %     figure;
%         %     for n = 1:N
%         %         subplot(R,R, n);
%         %         imagesc(shiftdim(clusteredData{ci}(n,:,:)).');
%         %         if N <= 9
%         %             colorbar;xlabel('Time');ylabel('Trials');title(['Neuron No. ' num2str(clustersinds{ci}(n))]);
%         %         end
%         % %         plot(shiftdim(clusteredData{ci}(n,:,:)))
%         %     end
%         disp(['Folder no. ' num2str(clusters(ci)) ': ' num2str(clustersinds{ci})]);
%     end
% end
eigsnum_row = 3;
eigsnum_col = 3;
eigsnum_trials = 3;
% [ err_rate, organized_data, row_order, col_order, trial_order ] = OrganizeData3D( data, data, row_aff, col_aff, trials_aff, row_perm, col_perm, trial_perm, eigsnum_col, eigsnum_row, eigsnum_trial );
row_thresh = 0.00;%0.03 for '8_17_14_1-45'  '8_17_14_46-80'

row_aff1 = threshold(row_aff, row_thresh);
[row_vecs, row_vals] = CalcEigs(row_aff1, eigsnum_row);
[ row_order ] = OrganizeDiffusion3DbyOneDim( data, row_vecs*row_vals );
[~, row_order]=sort(row_order);
PlotEmbedding(row_vecs*row_vals,row_order,'');
plotMeanOnFoldersAndTrials=true;

if plotMeanOnFoldersAndTrials
%     treeLevel = 1;
%     meanFolderTrials{treeLevel} = zeros(row_tree{treeLevel}.folder_count, nt);
%     for n = 1:length(row_order)
%         meanFolderTrials{treeLevel}(n, :) = mean(data(row_order(n), :, :), 3);
%     end
%     figure;imagesc(meanFolderTrials{treeLevel});colorbar;xlabel('Time');
    
    for treeLevel=1:length(row_tree)-1
        currData{treeLevel}=zeros(row_tree{treeLevel}.folder_count, nt, nT);
        currVar{treeLevel}=zeros(row_tree{treeLevel}.folder_count, nt);
        for tl = 1:row_tree{treeLevel}.folder_count
            currData{treeLevel}(tl, :, :) = mean(data(row_tree{treeLevel}.clustering == tl, :, :),1);
            if sum(row_tree{treeLevel}.clustering == tl) ~= 1            
            currVar{treeLevel}(tl, :) = mean(std(data(row_tree{treeLevel}.clustering == tl, :, :)),3);
            end
        end
        
        selinds = find(sum(currVar{treeLevel} < 0.3,2)==nt);
        
        curr_row_aff{treeLevel} = CalcEmdAff3D(permute(currData{treeLevel}, [2 3 1]), trials_tree, col_tree, params.col_emd, params.trials_emd, params.col_emd, ~params.init_aff_col.on_rows);
        [curr_row_vecs{treeLevel}, curr_row_vals{treeLevel}] = CalcEigs(curr_row_aff{treeLevel} , eigsnum_row);
        [ curr_row_order{treeLevel} ] = OrganizeDiffusion3DbyOneDim( currData{treeLevel}, curr_row_vecs{treeLevel}*curr_row_vals{treeLevel} );
        [~, curr_row_order{treeLevel}]=sort(curr_row_order{treeLevel});
        figure;
        PlotEmbedding(curr_row_vecs{treeLevel}*curr_row_vals{treeLevel},curr_row_order{treeLevel},'');
        meanFolderTrialsSel{treeLevel} = [];
        meanFolderTrials{treeLevel} = [];
        l=1;
        for n = 1:length(curr_row_order{treeLevel})
            meanFolderTrials{treeLevel}(n, :) = mean(currData{treeLevel}(curr_row_order{treeLevel}(n), :, :), 3);
            if all(currVar{treeLevel}(curr_row_order{treeLevel}(n), :) < 0.3)
            meanFolderTrialsSel{treeLevel}(l, :) = mean(currData{treeLevel}(curr_row_order{treeLevel}(n), :, :), 3);
            l = l + 1;            
            end
        end
        figure;subplot(2,1,1);imagesc(meanFolderTrials{treeLevel});colorbar;xlabel('Time');
        title(['Mean Neurons Response Tree Level: ' , num2str(treeLevel)]);
        subplot(2,1,2);imagesc(meanFolderTrialsSel{treeLevel});colorbar;xlabel('Time');
        title('Filtered By STD');
        if savefigs
%             print(fullfile(figspath1, ['FolderMeansLevel' num2str(clusterLevel)  '.pdf']),'-dpdf');
            saveas(gcf, fullfile(figspath1, ['FolderMeansLevel' num2str(treeLevel) '.jpg']),'jpg');
        end
    end
    
end


plotMeanOnTrials=false;
if plotMeanOnTrials
    
    meanMat = cell(length(row_tree), 1);
    for clusterLevel=length(row_tree):-1:2
        levelclustering = row_tree{clusterLevel}.clustering;
        clusters = unique(levelclustering);
        clusteredData = cell(length(clusters), 1);
        meanMat{clusterLevel} = cell(length(clusters),1);
        for ci = 1:length(clusters)
            
            clusteredData{ci} = data(levelclustering==clusters(ci), :, :);
            clustersinds{ci} = find(levelclustering==clusters(ci));
        end
        for ci = 1:length(clusters)
            orderingByEmbeding = row_order(clustersinds{ci});
            [~, IC] = sort(orderingByEmbeding);
            %         embd = row_vecs(clustersinds{ci},:)*row_vals;
            %         [ row_order1 ] = OrganizeDiffusion3DbyOneDim( clusteredData{ci}, embd );
            %         orderingByEmbeding = row_order1;
            %         [~, IC] = sort(orderingByEmbeding);
            N = sum(levelclustering==clusters(ci));
            if N==1
                figure;
                plot(mean(shiftdim(clusteredData{ci}(1,:,:)),2));
                xlabel('Time');
                title(['Tree Level ' num2str(clusterLevel) ' Cluster No. ' num2str(clusters(ci))]);
                set(gca,'YtickLabel',ylabelsstr) ;
                set(gca,'Ytick',sort(yloc, 'ascend')) ;
                
                axis tight;
            else
                figure;
                subplot(2,2,3);
                
                %             plotEmbedding(embd,IC,'');
                
                ism = ismember(1:size(row_vecs,1), clustersinds{ci});
                row_vecsInf = row_vecs;
                row_vecsInf(~ism,:) = inf;
                plotEmbedding(row_vecsInf*row_vals,row_order,'');
                
                subplot(2,1,1);
                l=1;ylabelsstr=cell(length(IC),1);
                yloc = fliplr(linspace(0,N/10, length(IC)));
                for n = IC
                    orderedN{clusterLevel}{ci}(l) = str2num(NeuronsLabels{clustersinds{ci}(n)});
                    meanMat{clusterLevel}{ci} = cat(2, meanMat{clusterLevel}{ci}, mean(shiftdim(clusteredData{ci}(n,:,:)),2));
                    %                 currLabels = NeuronsLabels(clustersinds{ci});
                    %                 currData = clusteredData{ci};
                    %                 orderedN{clusterLevel}{ci}(l) = str2double(currLabels{n});
                    %                 meanMat{clusterLevel}{ci} = cat(2, meanMat{clusterLevel}{ci}, mean(shiftdim(currData(n,:,:)),2));
                    
                    plot(yloc(l) + mean(shiftdim(clusteredData{ci}(n,:,:)),2));
                    yloc(l) = yloc(l) + mean(mean(shiftdim(clusteredData{ci}(n,:,:)),2));
                    %                 plot(yloc(l) + mean(shiftdim(currData(n,:,:)),2));
                    %                 yloc(l) = yloc(l) + mean(mean(shiftdim(currData(n,:,:)),2));
                    hold on;
                    ylabelsstr{l} = NeuronsLabels{clustersinds{ci}(n)};
                    
                    %                 ylabelsstr{l} = currLabels{n};
                    l=l+1;
                    %         %         subplot(R,R, n);
                    %         %         imagesc(shiftdim(clusteredData{ci}(n,:,:)).');title(['ROI No. ' NeuronsLabels{clustersinds{ci}(n)}]);
                    %         %         if N <= 9
                    %         %             colorbar;xlabel('Time');ylabel('Trials');
                    %         %         end
                    %
                    %         %                 plot(shiftdim(clusteredData{ci}(n,:,:)))
                    %         %         errorbar(1:155, mean(shiftdim(clusteredData{ci}(n,:,:))),std(shiftdim(clusteredData{ci}(n,:,:))))
                end
                xlabel('Time');
                
                set(gca,'YtickLabel',ylabelsstr) ;
                set(gca,'Ytick',sort(yloc, 'ascend')) ;
                
                axis tight;
                title(['Tree Level ' num2str(clusterLevel) ' Cluster No. ' num2str(clusters(ci))]);
                subplot(2,2,4);
                imagesc(meanMat{clusterLevel}{ci}.');
                set(gca,'YtickLabel',ylabelsstr) ;
                set(gca,'Ytick',1:length(yloc)) ;
                ylabel('Neurons');
                xlabel('Time');
                title(['Tree Level ' num2str(clusterLevel) ' Cluster No. ' num2str(clusters(ci))]);
            end
            if savefigs
                print(fullfile(figspath1, ['Level' num2str(clusterLevel) 'Cluster' num2str(clusters(ci)) '.pdf']),'-dpdf');
                saveas(gcf, fullfile(figspath1, ['Level' num2str(clusterLevel) 'Cluster' num2str(clusters(ci)) '.jpg']),'jpg');
            end
            
        end
        %     ca;
    end
end
% save(wrkspname);


