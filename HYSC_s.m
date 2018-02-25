function output = HYSC_s(data, varargin)

% run Number of sIngle-cell Clusters Estimation (HYSC) to automatically
% estimate the number of cells clusters. 
% Input:
% data - M by N matrix, M is the number of variables (genes) and N is the number of cells.
% varargin  -  'ClusterNum' - the given number of clusters
%              'Nmax'       - the maximum number of cell clusters
%              'minClusterSize_number' - the minimum number of samples (cells) in each cluster
%              'dimClustering' - the number of principle components adopted by k-means clustering
%              nc - the number of correlated component in X of OPLSDA model
%              ncox - the number of orthogonal component in X of OPLSDA model
%              ncoy - the number of correlated component in Y of OPLSDA model
%              r2cutoff - OPLSDA parameters, r2 cutoff value to identify discriminatory variables
%              pcutoff - OPLSDA parameters, p cutoff value to identify discriminatory variables
%              preprocessing - OPLSDA parameters, pre processing methods, 'uv', 'pa', ''mc', by default 'uv'.
%
% Output:
% the structure 'output'
% 'output.groupIdx' - the cells categorization result based on the estimated number of cells clusters
% 'output.biomarkers' - the variable genes identified by hybrid clustering
%
% Author Xin Zou & Jie Hao, SJTU, China
% Copyright Xin Zou & Jie Hao



varList = {'ClusterNum', 'Nmax', 'minClusterSize_number', 'dimClustering','nc', 'ncox', 'ncoy', 'r2cutoff', 'pcutoff', 'preprocessing', 'perp'};
convgThrh = {0, 10, 10, 5, 1, 1, 0, 0.5, 0.05, 'uv', 30}; % default values
if ~isempty(varargin)
    L = length(varargin);
    if L == 1
        varargin = varargin{1};
        L = length(varargin);
    end
    for i = 1:2:L-1
        idx = find(strcmp(varargin{i}, varList));
        convgThrh{idx} = varargin{i+1};
    end
end

ClusterNum = convgThrh{1};
Nmax = convgThrh{2};
minClusterSize_number = convgThrh{3};
dim = convgThrh{4};
oplsda_para.nc = convgThrh{5};
oplsda_para.ncox = convgThrh{6};
oplsda_para.ncoy = convgThrh{7};
oplsda_para.r2cutoff = convgThrh{8};
oplsda_para.pcutoff = convgThrh{9};
oplsda_para.preprocessing = convgThrh{10};
perp = convgThrh{11};
oplsda_para.NP = 0;

%% preprocessing
temp = sum(data~=0,2);
Idx = find(temp>1);
data = data(Idx,:); % remove the genes with all zeros reads?
var = 1:size(data,1);
data = data';
data(data==0) = 10^(-6)*rand(1, length(find(data==0))); % add a small perturbation
%% HYSC
N = 2; % initial number of clusters
biomarkers_idx_pool_ori2 = [];
biomarkers_idx_pool_unsigned = 1:size(data,2);
randVar = [];
% % negative control to discriminate the cases of N=1 from N>1
while 1 % identify variable genes
    idx_biomarkers = unique(biomarkers_idx_pool_unsigned);
    % % PCA & k-means
    X = data(:,[idx_biomarkers, randVar]);
    X_temp = (X-ones(size(X, 1),1)*mean(X,1));
    [~,SCORE] = pca(X_temp);
    temp = SCORE(:,1:min(dim,size(SCORE,2)));
    IDX = kmeans(temp,N, 'replicates', 20);
    
    % pairwise comparison between clusters
    biomarkers_idx_pool_unsigned = [];
    C = combnk(1:N,2);
    P_temp = ones(size(C,1)+N, size(data,2));
    randVar = [];
    for i = 1:size(C,1) % for pair of clusters
        tempData1 = data(IDX==C(i,1), :);
        tempData2 = data(IDX==C(i,2), :);
        model = OPLSDA(tempData1, tempData2, var, oplsda_para);
        [~, randTemp] = sort(model.W2);
        randTemp = randTemp(1:floor(length(model.sig_idx)*0.05)+1);
        if isempty(model) || isempty(model.sig_idx) || length(model.sig_idx)<2 || length(find(IDX==C(i,1)))<minClusterSize_number || length(find(IDX==C(i,2)))<minClusterSize_number
            output.groupIdx = ones(size(IDX));
            return;
        end
        
        biomarkers_idx_pool_unsigned = [biomarkers_idx_pool_unsigned, model.sig_idx];
        randVar = [randVar, randTemp];
        P_temp(i,:) = model.P;
    end
    
    
    % convergence detection
    biomarkers_idx_pool_unsigned = unique(biomarkers_idx_pool_unsigned);
    randVar = unique(randVar);
    if ~isempty(biomarkers_idx_pool_ori2)
        c = intersect(biomarkers_idx_pool_ori2, biomarkers_idx_pool_unsigned);
        if (length(biomarkers_idx_pool_unsigned)-length(c))/length(biomarkers_idx_pool_ori2)<0.05 % if the variable genes identified in two consecutive iterations have more than 95% overlap, converge
            break;
        end
    end
    biomarkers_idx_pool_ori2 = biomarkers_idx_pool_unsigned;
end

% negative control, only applied for the cluster <5000 cells
if size(X,1)<5000
    X = data(:,biomarkers_idx_pool_unsigned);
    [COEFF,SCORE,latent] = pca(X);
    ydata = tsne(X,[],SCORE(:,1:min(dim, length(model.sig_idx))), perp);
    [~, p] = ttest2(ydata(IDX==1, 1), ydata(IDX==2, 1));
    if p>=0.05
        output.groupIdx = ones(size(IDX));
        return;
    end
end

if N == min(Nmax, ClusterNum) % if the maximum number of clusters reached, the algorithm terminates
    groupIdx = IDX;
    output.groupIdx = groupIdx;
    output.biomarkers = biomarkers_idx_pool_global;
    return;
end
biomarkers_idx_pool_ori = [];
biomarkers_idx_pool_global = unique(biomarkers_idx_pool_unsigned);
N = N+1;

% % Hybrid-clustering
while 1 % loop for estimating N
    if N>2 % save the resutls obtained in the previous iteration
        IDX_lastIte = IDX;
        biomarkers_idx_pool_ori = biomarkers_idx_pool_global;
    end
    
    biomarkers_idx_pool_ori2 = [];
    while 1 % loop for variable genes identification
        idx_biomarkers = unique(biomarkers_idx_pool_unsigned);
        % % PCA & k-means
        X = data(:,[idx_biomarkers, randVar]);
        X_temp = (X-ones(size(X, 1),1)*mean(X,1));
        [~,SCORE] = pca(X_temp);
        temp = SCORE(:,1:min(dim,size(SCORE,2)));
        IDX = kmeans(temp,N, 'replicates', 20);
        
        a = unique(IDX);
        out = [a,histc(IDX(:),a)]; % the number of samples in each cell cluster
        if min(out(:,2))<minClusterSize_number && ClusterNum ==0 % if the size of a cluster is too small, converges.
            groupIdx = IDX_lastIte;
            output.groupIdx = groupIdx;
            output.biomarkers = biomarkers_idx_pool_global;
            return;
        end
        
        % pairwise comparison between clusters
        biomarkers_idx_pool_unsigned = [];
        C = combnk(1:N,2);
        P_temp = ones(size(C,1)+N, size(data,2));
        biomarkers_cell = cell(1, size(C,1));
        data_cluster_idx_cell = cell(size(C));
        randVar = [];
        for i = 1:size(C,1) % compare each pair of cell clusters
            tempData1 = data(IDX==C(i,1), :);
            tempData2 = data(IDX==C(i,2), :);
            model = OPLSDA(tempData1, tempData2, var, oplsda_para);
            [~, randTemp] = sort(model.W2);
            randTemp = randTemp(1:floor(length(model.sig_idx)*0.05)+1);
            if isempty(model) % if the two clusters are unseparable, converge
                groupIdx = IDX_lastIte;
                biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
                output.groupIdx = groupIdx;
                output.biomarkers = biomarkers_idx_pool_global;
                return;
            end
            biomarkers_idx_pool_unsigned = [biomarkers_idx_pool_unsigned, model.sig_idx];
            randVar = [randVar, randTemp];
            P_temp(i,:) = model.P;
            biomarkers_cell{i} = model.sig_idx;
            data_cluster_idx_cell{i, 1} = find(IDX==C(i,1));
            data_cluster_idx_cell{i, 2} = find(IDX==C(i,2));
        end
        for i = 1:N % % compare each cell cluster to the rest of cells
            tempData1 = data(IDX==i, :);
            tempData2 = data(IDX~=i, :);
            model = OPLSDA(tempData1, tempData2, var, oplsda_para);
            [~, randTemp] = sort(model.W2);
            randTemp = randTemp(1:floor(length(model.sig_idx)*0.05)+1);
            if isempty(model)% if the two clusters are unseparable, converge
                groupIdx = IDX_lastIte;
                biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
                output.groupIdx = groupIdx;
                output.biomarkers = biomarkers_idx_pool_global;
                return;
            end
            P_temp(size(C,1)+1,:) = model.P;
            biomarkers_idx_pool_unsigned = [biomarkers_idx_pool_unsigned, model.sig_idx];
            randVar = [randVar, randTemp];
        end
        P_min = min(P_temp);
        
        % variable genes identification iteration termination
        biomarkers_idx_pool_unsigned = unique(biomarkers_idx_pool_unsigned);
        if ~isempty(biomarkers_idx_pool_ori2)
            c = intersect(biomarkers_idx_pool_ori2, biomarkers_idx_pool_unsigned);
            if (length(biomarkers_idx_pool_unsigned)-length(c))/length(biomarkers_idx_pool_ori2)<0.1 % if the variable transcripts identified in two consecutive iterations have more than 90% overlap, converge
                break;
            end
        end
        biomarkers_idx_pool_ori2 = biomarkers_idx_pool_unsigned;
    end
    
    % N estimation iteration termination
    if ~isempty(find(cellfun('isempty',biomarkers_cell))) && ClusterNum ==0 % if no variable gene observed in any comparison, converge
        groupIdx = IDX_lastIte;
        biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
        output.groupIdx = groupIdx;
        output.biomarkers = biomarkers_idx_pool_global;
        return;
    end
    if N == Nmax  && ClusterNum ==0 || N == ClusterNum % if the maximum number of clusters reached or N has reached the given number of clusters, converge
        groupIdx = IDX;
        biomarkers_idx_pool_global = unique(biomarkers_idx_pool_global);
        output.groupIdx = groupIdx;
        output.biomarkers = biomarkers_idx_pool_global;
        return;
    end
    
    % % gene expressin profile binaryzatin
    biomarkers_idx_pool_global = unique([biomarkers_idx_pool_unsigned, biomarkers_idx_pool_global]);
    SIGMA = binaryzation(data, IDX, biomarkers_idx_pool_global);
    % if there exist dulplicate columns in SIGMA or the matrix 'SIGMA' is empty, converge
    if ~isempty(SIGMA)
        D = pdist(SIGMA', 'hamming');
        if min(D) == 0
            groupIdx = IDX_lastIte;
            biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
            output.groupIdx = groupIdx;
            output.biomarkers = biomarkers_idx_pool_global;
            return;
        end
    else
        groupIdx = IDX_lastIte;
        biomarkers_idx_pool_global = unique(biomarkers_idx_pool_ori);
        output.groupIdx = groupIdx;
        output.biomarkers = biomarkers_idx_pool_global;
        return;
    end
    
    
    N=N+1; % update the putative number of cell clusters
    
end
