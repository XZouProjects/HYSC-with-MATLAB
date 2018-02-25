function output = HYSC()

% run Hybrid-clustering for Single-cell Categrization (HYSC)
%
% Input parameters:
% dataPath                 - The directory of .xlsx file containing expression data. The first column contains the gene symbles and the first row contains cell IDs. The dropouts should be zero-inflated.
% ClusterNum               - the given number of clusters. If ClusterNum = 0, the algorithm automatically estimates the number of clusters
% Nmax                     - the maximum number of cell clusters. By default, Nmax = 30
% minClusterSize_number    - the minimum number of samples (cells) in each cluster. By default, minClusterSize_number = 10
% dimClustering            - the number of principle components adopted by k-means clustering. By default, dimClustering = 30
% r2cutoff                 - OPLSDA parameters, r2 cutoff value to identify discriminatory variables. By default 0.5
% pcutoff                  - OPLSDA parameters, p cutoff value to identify discriminatory variables. By default 0.05 after bonferroni correction
% perp                     - perplexity parameter of tSNE. By default, perp = 30
% maxHYSCLayer             - the maximum number of HYSC layers. By default, maxHYSCLayer = 5
% cores                    - the number of cores used in parallel computation. By default, cores = 8
% tSNEScores               - the index of tSNE scores adopted for scores plot. By default, the first two components are illustrated.
% geneFiltering_cellCounts - the mimimum number of cells in which a gene must express. By default >0.
% geneFiltering_var        - the mimimum variance of a valid gene. By default >0.
%
% Output:
% The output is stored in an excel file with 3 sheets, 
% Sheet 1 - the cell categorization indexs
% Sheet 2 - the list of variable genes
% Sheet 3 - the tSNE scores, cells vs tSNE components
%
% Author Xin Zou & Jie Hao, SJTU, China
% Copyright Xin Zou & Jie Hao

warning off;

varList = {'ClusterNum', 'Nmax', 'minClusterSize_number', ...
    'dimClustering', 'r2cutoff', 'pcutoff', 'perp', 'maxHYSCLayer', 'cores', 'tSNEScores', 'geneFiltering_cellCounts', 'geneFiltering_var'};
convgThrh = {0, 30, 10, 5, 0.5, 0.05, 30, 5, 8, 1:2, 0, 0}; % default values


fid = fopen('HYSC Input.txt', 'r');
Inputs = [];
i = 1;
while ~feof(fid)
    S=fgetl(fid);
    tf = find(isspace(S));
    if i == 1
        dataPath = S(tf(end)+1:end);
    else
        Inputs{i-2} = S(1:tf(end)-2);
        Inputs{i-1} = str2double(S(tf(end)+1:end));
    end
    i = i+2;
end
fclose(fid);

if ~isempty(Inputs)
    L = length(Inputs);
    for i = 1:2:L-1
        idx = find(strcmp(Inputs{i}, varList));
        convgThrh{idx} = Inputs{i+1};
    end
end

ClusterNum = convgThrh{1}; % the number of cell clusters if specified
Nmax = convgThrh{2}; % the maximum number of cell clusters
minClusterSize_number = convgThrh{3}; % the minimum number of samples (cells) in each cluster
dim = convgThrh{4}; % the number of principle components adopted by k-means clustering
oplsda_para.r2cutoff = convgThrh{5}; % the OPLSDA loading coefficient cutoff for variable genes identification
oplsda_para.pcutoff = convgThrh{6}; % the OPLSDA p-value cutoff for variable genes identification
perp = convgThrh{7}; % perplexity parameter of tSNE
maxHYSCLayer = convgThrh{8}; % the maximum number of HYSC layers
cores = convgThrh{9}; % the number of cores used in parallel computation
tSNEScores = convgThrh{10}; % the index of tSNE scores adopted for scores plot
geneFiltering_cellCounts = convgThrh{11}; % 
geneFiltering_var = convgThrh{12}; % 
oplsda_para.nc = 1; % the number of correlated component in X of OPLSDA model
oplsda_para.ncox = 1; % the number of orthogonal component in X of OPLSDA model
oplsda_para.ncoy = 0; % the number of correlated component in Y of OPLSDA model
oplsda_para.preprocessing = 'uv'; % the scaling method in OPLSDA is 'unit variance'
oplsda_para.NP = 0; % no permutation test is performed

%% Parameters setup
paraList_HYSC = {'ClusterNum',ClusterNum, 'Nmax',Nmax, 'minClusterSize_number', minClusterSize_number,...
    'dimClustering', dim, 'nc', oplsda_para.nc, 'ncox', oplsda_para.ncox, ...
    'ncoy', oplsda_para.ncoy, 'r2cutoff',oplsda_para.r2cutoff...
    'pcutoff', oplsda_para.pcutoff, 'preprocessing', oplsda_para.preprocessing, 'perp', perp};


%% load data
[~, ~, raw] = xlsread(dataPath);
data = cell2mat(raw(2:end, 2:end));
geneList = raw(2:end, 1);

Idx = [];
for i = 1:size(data,1)
    num = length(find(data(i,:)~=0));
    if num>geneFiltering_cellCounts
        if var(data(i,:))>geneFiltering_var
            Idx = [Idx, i];
        end
    end
end
data = data(Idx,:);
geneList = geneList(Idx);



% HYSC cluster symbles
types = cell(1,26);
s=(1:26)+96;
str=char(s);
for i = 1:length(str)
    types{i} = str(i);
end

clusterMarkParal = cell(1,cores);
unsepClusterMarkParal = cell(1,cores);
%% parallel computing
startmatlabpool(cores);
biomarkers = cell(1, cores);
parfor paral = 1:cores
% for paral = 1:cores
    clusterMark = cell(size(data, 2), 1);
    clusterMark(:) = cellstr(types{1});
    unsepClusterMark = [];
    while 1
        uniqueClusterMark = unique(clusterMark);
        NumOfClusters = length(uniqueClusterMark);
        L = zeros(1,NumOfClusters);
        for i = 1:NumOfClusters
            L(i) = length(char(uniqueClusterMark{i}));
        end
        currHYSCLayer = max(L); % the current HYSC layer
        if NumOfClusters < Nmax && currHYSCLayer <= maxHYSCLayer
            for n = 1:NumOfClusters % further categorize each cell cluster
                if ismember(uniqueClusterMark{n}, unsepClusterMark) == 0 % check if the subset is separable
                    idx = find(ismember(clusterMark, uniqueClusterMark{n}));
                    if length(idx) >= minClusterSize_number % the size of each cell cluster should be no smaller than the minimum cluster size
                        dataBuffer = data(:, idx);
                        temp = sum(dataBuffer~=0,2);
                        Idx = find(temp>2);
                        dataBuffer = dataBuffer(Idx,:);
                        outputHYSC = HYSC_s(dataBuffer, paraList_HYSC); % categorize the data subset using HYSC
                        if max(outputHYSC.groupIdx) == 1 % if the subset is unseparable
                            unsepClusterMark = [unsepClusterMark, uniqueClusterMark(n)];
                        else % assign a label to each obtained cell cluster
                            for j = 1:max(outputHYSC.groupIdx) 
                                idx2 = find(outputHYSC.groupIdx==j);
                                for k = 1:length(idx2)
                                    clusterMark{idx(idx2(k))} = [clusterMark{idx(idx2(k))}, types{j}];
                                end
                            end
                            biomarkers{paral} = [biomarkers{paral}; Idx(outputHYSC.biomarkers)];
                        end
                    else
                        unsepClusterMark = [unsepClusterMark, uniqueClusterMark(n)];
                    end
                end
            end
        end
        if length(unsepClusterMark) == NumOfClusters && NumOfClusters == 1
            break;
        end
        if ~isempty(unsepClusterMark)  && length(unsepClusterMark) == NumOfClusters || NumOfClusters >= Nmax || currHYSCLayer == maxHYSCLayer
            break;
        end
        
    end
    unsepClusterMarkParal{paral} = unsepClusterMark;
    clusterMarkParal{paral} = clusterMark;
    biomarkers{paral} = unique(biomarkers{paral});
end

% identify the globle optima
NumOfClusters_paral = zeros(1,cores);
clusterIdx = cell(1,cores);
uniqueClusterMark = cell(1,cores);
for paral = 1:cores
    uniqueClusterMark{paral} = unique(clusterMarkParal{paral});
    NumOfClusters_paral(paral) = length(uniqueClusterMark{paral});
    clusterIdx{paral} = cell(1,NumOfClusters_paral(paral));
    for j = 1:NumOfClusters_paral(paral)
        idx = find(strcmp(uniqueClusterMark{paral}{j},clusterMarkParal{paral}));
        clusterIdx{paral}{j} = idx;
    end
end

[freqs, centers] = hist(NumOfClusters_paral, unique(NumOfClusters_paral));
[~, IX] = max(freqs);
idx = find(NumOfClusters_paral==centers(IX));
idx = idx(1);


clusterMark = clusterMarkParal{idx};
unsepClusterMark = unsepClusterMarkParal{idx};
NumOfClusters = NumOfClusters_paral(idx);
biomarkers_final = biomarkers{idx};
closematlabpool;

%% the output results

if length(unsepClusterMark) == NumOfClusters && NumOfClusters == 1 % if the dataset is unseparable
    output.biomarkers = [];
    output.clusterIdx = ones(size(data, 2), 1);
    currentFolder = pwd;
    delete([currentFolder, '\HYSC Output.xlsx']);
    xlswrite([currentFolder, '\HYSC Output.xlsx'], output.clusterIdx);
else
    uniqueClusterMark = unique(clusterMark);
    NumOfClusters = length(uniqueClusterMark);
    IDX = zeros(size(clusterMark));
    
    for n = 1:NumOfClusters
        IDX(ismember(clusterMark, uniqueClusterMark{n})) = n;
    end
    
    output.clusterIdx = IDX;
    output.biomarkers = geneList(biomarkers_final);
    currentFolder = pwd;
    delete([currentFolder, '\HYSC Output.xlsx']);
    xlswrite([currentFolder, '\HYSC Output.xlsx'], output.clusterIdx, 1);
    xlswrite([currentFolder, '\HYSC Output.xlsx'], output.biomarkers, 2);
%% illustrate the clustering patterns in tSNE scatter plot
% % only the clustering results obtained in the up to first 3 HYSC layers are illustrated
    X = data(biomarkers_final,:)';
    [COEFF,SCORE,latent] = pca(X);
    if size(X,1)<5000 % for the dataset smaller than 5000 cells, used tSNE scores plot
        ydata = tsne(X,[],SCORE(:,1:min(5,size(X,2))), perp);
    else % otherwise pca scores plot
        ydata = SCORE(:,1:min(5,size(X,2)));
    end
    output.tSNE = ydata;
    xlswrite([currentFolder, '\HYSC Output.xlsx'], output.tSNE, 3);
    e = actxserver('Excel.Application'); % # open Activex server
    ewb = e.Workbooks.Open([currentFolder, '\HYSC Output.xlsx']); % # open file (enter full path!)
    ewb.Worksheets.Item(1).Name = 'Cluster ID'; % # rename 1st sheet
    ewb.Worksheets.Item(2).Name = 'Variable Genes'; % # rename 1st sheet
    ewb.Worksheets.Item(3).Name = 'tSNE scores'; % # rename 1st sheet
    ewb.Save % # save to the same file
    ewb.Close(false)
    e.Quit

    
    figure;
    subplot('Position',[0.14 0.16 0.73 0.76]);
    colormap(jet);
    Color=get(gcf,'Colormap');
    
    % color code the cells according to HYSC hierarchy
    layerInfor = cellfun('length',clusterMark)-1;
    cate_layerd = zeros(length(layerInfor), max(layerInfor));
    for j = 1:max(layerInfor) % for each layer
        buffer = cell(1,length(layerInfor));
        for jj = 1:length(layerInfor)
            if layerInfor(jj)>=j % only apply for the cluters obtained by current layer or subsequent layers
                buffer{jj} = clusterMark{jj}(j+1); % exclude the firt psedo layer category
            end
        end
        idx1 = find(~cellfun('isempty', buffer)); % identify the index of none empty 'buffer' element
        buffer_unique = unique(buffer(idx1));
        temp_idx = cell(length(buffer_unique));
        for jj = 1:length(buffer_unique)
            idx2 = find(strcmp(buffer_unique{jj}, buffer(idx1)));
            cate_layerd(idx1(idx2),j) = jj;
        end
    end
    % assign a color to each cell
    color_cate_layerd = cell(size(cate_layerd,1), min(max(layerInfor),3));
    for j = 1:min(max(layerInfor),3) % for each layer
        if j == 1 % the 1st layer
            buffer = unique(cate_layerd(:,j));
            linIdx = floor(linspace(1,64,length(buffer)+1)); % linearly spaced index vector
            for jj = 1:length(buffer)
                color_cate_layerd(find(cate_layerd(:,j)==buffer(jj)),j) = {[linIdx(jj),linIdx(jj+1)]};
            end
        elseif j == 2 % the 2nd layer
            buffer = unique(cate_layerd(:,1)); % for the 1st layer
            for jj = 1:length(buffer)
                idx = find(cate_layerd(:,j-1)==buffer(jj)); % check if the upper layer has subsequent layers
                if sum(cate_layerd(idx,j))>0 % there is a subsequent layer
                    buffer2 = unique(cate_layerd(idx,j));
                    linIdx = floor(linspace(color_cate_layerd{idx(1), j-1}(1), color_cate_layerd{idx(1), j-1}(2), length(buffer2)+1)); % linearly spaced index vector according the color range of the upper layer
                    %                     linIdx_cell = cell(1, length(buffer2));
                    for jjj = 1:length(buffer2)
                        idx2 = find(cate_layerd(idx,j)==buffer2(jjj));
                        color_cate_layerd(idx(idx2),j) = {[linIdx(jjj),linIdx(jjj+1)]};
                    end
                end
            end
        elseif j == 3 % the 3rd layer
            buffer1 = unique(cate_layerd(:,1)); % the 1st layer
            buffer2 = unique(cate_layerd(:,2)); % the 2nd layer 
            for jj = 1:length(buffer1)
                for jjj = 1:length(buffer2)
                    idx = find(cate_layerd(:,1)==buffer1(jj) & cate_layerd(:,2)==buffer2(jjj)); % check if the upper layer has subsequent layers
                    
                    if sum(cate_layerd(idx,j))>0 % there is a subsequent layer
                        buffer3 = unique(cate_layerd(idx,j));
                        linIdx = floor(linspace(color_cate_layerd{idx(1), j-1}(1), color_cate_layerd{idx(1), j-1}(2), length(buffer3)+1)); % linearly spaced index vector according the color range of the upper layer
                        %                     linIdx_cell = cell(1, length(buffer2));
                        for j4 = 1:length(buffer3)
                            idx2 = find(cate_layerd(idx,j)==buffer3(j4));
                            color_cate_layerd(idx(idx2),j) = {[linIdx(j4),linIdx(j4+1)]};
                        end
                    end
                    
                end
            end
        end
        
    end
    
    for i = 1:length(IDX)
        idx = max(find(~cellfun('isempty', color_cate_layerd(i,:))));
        color_code(i) = round(mean(color_cate_layerd{i,idx}));
        plot(ydata(i,tSNEScores(1)), ydata(i,tSNEScores(2)), '.', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 0.5);
        t = text(ydata(i,tSNEScores(1)), ydata(i,tSNEScores(2)),num2str(IDX(i)));
        t.Color = Color(color_code(i),:);
        t.FontSize = 15;
        if i == 1
            hold on;
        end
    end
    set(gca, 'fontsize', 14);
    xlabel('tSNE_1');
    ylabel('tSNE_2');
    
    subplot('Position',[0.9 0.16 0.05 0.76]);
    imagesc(1,1:-1/64:0,(1:64)');
    set(gca,'xTickLabel','');
    set(gca,'yTickLabel','');
    
end






























