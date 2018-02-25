function SIGMA = binaryzation(data, groupIdx, biomarkersIdx)

% run gene expression profiles binaryzation
% Inputs:
%              'data'- the data matrix with genes being columns and cells being rows
%              'groupIdx' - the cell categorization results obtained by HYSIC
%              'biomarkersIdx' - the gene markers index obtained by HYSIC
% Outputs:
%              'SIGMA' - the binaryzed expression matrix
%
% Author Xin Zou & Jie Hao, SJTU, China
% Copyright Xin Zou & Jie Hao

% binarize the expression profiles of each variable gene
X_binary = zeros(length(groupIdx), length(biomarkersIdx));
for i = 1:length(biomarkersIdx)
    temp = data(:,biomarkersIdx(i));
    IDX = kmeans(temp, 2, 'replicates', 20);
    mu(1) = mean(temp(IDX == 1));
    mu(2) = mean(temp(IDX == 2));
    [~, IX] = max(mu);
    thresh = min(temp(IDX == IX));
    X_binary(temp>=thresh,i) = 1;
    X_binary(temp<thresh,i) = 0;
end
clear P;

% Construct the binary matrix SIGMA
for i = floor(log2(max(groupIdx))):2^max(groupIdx) 
    if i == 1 % when there is only one gene group
        % contingency table
        counts = zeros(2,max(groupIdx));
        for j = 1:max(groupIdx) % for each cell category
            counts(1,j) = length(find(X_binary(groupIdx==j,:)==1));
            counts(2,j) = length(find(X_binary(groupIdx==j,:)==0));
%             grouped_genes{1,j} = X_binary(groupIdx==j,:)';
        end
        % hypergeometric test
        K = sum(counts);
        n = sum(counts,2);
        M = sum(K);
        for k = 1:size(counts,1)
            P{i}(k,:) = 1-cdf('hyge',counts(k,:), M, K, n(k));
        end
        % construct binary metrix SIGMA
        temp = P{i}(1,:);
        temp(find(P{i}(1,:)<0.05)) = 1;
        temp(find(P{i}(1,:)>=0.05)) = 0;
    else
        IDX = kmeans(X_binary', i, 'replicates', 20, 'Distance', 'hamming');
        temp = [];
        geneNgroup = zeros(1,i);
        for jj = 1:i % for each gene group
            % contingency table
            counts = zeros(2,max(groupIdx));
            for j = 1:max(groupIdx) % for each cell category
                counts(1,j) = length(find(X_binary(groupIdx==j,IDX==jj)==1));
                counts(2,j) = length(find(X_binary(groupIdx==j,IDX==jj)==0));
%                 grouped_genes{jj,j} = X_binary(groupIdx==j,IDX==jj)';
            end

            % hypergeometric test
            K = sum(counts);
            n = sum(counts,2);
            M = sum(K);
            for k = 1:size(counts,1)
                P{i,jj}(k,:) = 1-cdf('hyge',counts(k,:), M, K, n(k));
            end
            % construct binary metrix SIGMA
            temp_P = P{i,jj}(1,:);
            temp_P(find(P{i,jj}(1,:)<0.05)) = 1;
            temp_P(find(P{i,jj}(1,:)>=0.05)) = 0;
            temp = [temp; temp_P];
            geneNgroup(jj) = length(find(IDX==jj));
        end
        % convergence detection
        D = pdist(temp, 'hamming');
        if min(D) == 0 || min(geneNgroup)<10
            if exist('temp_ori')
                SIGMA = temp_ori; % rows are gene clusters and columns are cells clusters
            else
                SIGMA = [];
            end
            return;
        end
    end
    temp_ori = temp;
end