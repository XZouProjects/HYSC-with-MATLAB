function model = OPLSDA(data1,data2, geneIdx, oplsda_para)

% run OPLS-DA and display results
% Input:
% "data1" and "data2" are groups of spectra from two biological classes
% "geneIdx" is the indexs of variables (genes)

% Output:
% "model" contains all OPLS-DA results, where the discriminatory variables are in the model.sig_idx.

% The parameters are set in the OPLSDA para.txt file
% nrcv - number of folds cross-validation
% nc - number of correlated variables in X
% ncox - number of orthogonal variables in X
% ncoy - number of orthogonal variables in Y
% r2cutoff - r2 cutoff value to identify discriminatory variables
% pcutoff - p cutoff value to identify discriminatory variables
% np - number of permutation tests
% preprocessing - pre processing methods, 'uv', 'pa', ''mc'

nc = oplsda_para.nc;
ncox = oplsda_para.ncox;
ncoy = oplsda_para.ncoy;
r2cutoff = oplsda_para.r2cutoff;
pCutoff = oplsda_para.pcutoff;
prep = oplsda_para.preprocessing;
pCutoff = pCutoff/length(geneIdx);
NP = oplsda_para.NP;

%% format data
X = [data1;data2];
Y = zeros(size(X,1),2);
Y(1:size(data1),1) = 1;
Y(size(data1)+1:end,2) = 1;

%%
if NP == 0
    model=o2pls_m(X,Y,prep,nc,ncox,ncoy,0, r2cutoff, pCutoff);
    if isempty(model)
        return;
    end
end