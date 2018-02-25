function model=o2pls_m(X,Y,prep,nc,ncox,ncoy,nrcv, r2cutoff, pCutoff)
%
% O2-PLS function used to build the model
% For further predictions use the o2pls_pred function
%
% Input:
% X - X matrix input
% Y - Y matrix input
% prep - preprocessing available
%           - 'no' : no preprocessing
%           - 'mc' : meancentering
%           - 'uv' : Univariance scaling
%
% nc - number of correlated components, it can be determinated by PCA of Y'X
% ncox - number of orthogonal components in X
% ncoy - number of orthogonal component in Y
% nrcv - number of fold in the cross-validation (full cross validation)
% v1.0, Xin Zou, Kent

opt = 1;
BT = 0;

if nargin<8
    r2cutoff=0.3;
    pCutoff=0.05;
end


muClass1 = mean(X((Y(:,1)==1),:));
muClass2 = mean(X((Y(:,2)==1),:));
Signs = sign(muClass2-muClass1);

if nargin < 7
    nrcv = 0; % no cross validation
end

[nsx,nvx]=size(X);
[nsy,nvy]=size(Y);

model.ns = nsx;

model.nc=nc;
model.ncox=ncox;
model.ncoy=ncoy;

if nsx~=nsy
    disp('Number of samples are different in X and Y')
    model=[];
    return
end

model.MeanX=mean(X);
model.MeanY=mean(Y);

model.StdX=std(X);
model.StdY=std(Y);

model.SSX=sum(sum(X.^2));
model.SSY=sum(sum(Y.^2));
model.CSX=sum(X.^2);
model.CSY=sum(Y.^2);
% Preprocessing

model.modelType = 'pls';



if nrcv == 0 % no CV
    switch prep
        
        case 'no'
            model.preprocessing='no';
            
        case 'mc'
            %         disp('Meancentering')
            model.preprocessing='mc';
            X=X-repmat(model.MeanX,nsx,1);
            Y=Y-repmat(model.MeanY,nsy,1);
            
        case 'uv'
            % disp('Univariance scaling')
            model.preprocessing='uv';
            X=X-repmat(model.MeanX,nsx,1);
            Y=Y-repmat(model.MeanY,nsy,1);
            X=X./repmat(model.StdX,nsx,1);
            Y=Y./repmat(model.StdY,nsy,1);
            
        case 'pa'
            model.preprocessing='pa';
            X=X-repmat(model.MeanX,nsx,1);
            Y=Y-repmat(model.MeanY,nsy,1);
            X=X./repmat(sqrt(model.StdX),nsx,1);
            Y=Y./repmat(sqrt(model.StdY),nsy,1);
            
        otherwise
            disp('Unknown Preprocessing')
            model=[];
            return
            
    end
    X(find(isnan(X))) = 0;
    M=O2pls(X,Y,nc,ncox,ncoy,'standard');
    if isempty(M)
        return;
    end
    M.preprocessing = model.preprocessing;
    %     M.W = M.W;
    M.ns = model.ns;
    
    M.nc = model.nc;
    M.ncox = model.ncox;
    M.ncoy = model.ncoy;
    
    M.MeanX = model.MeanX;
    M.MeanY = model.MeanY;
    
    M.StdX = model.StdX;
    M.StdY = model.StdY;
    M.X_out = X;
    
    M.SSX = model.SSX;
    M.SSY = model.SSY;
    
    M.CSX = model.CSX;
    M.CSY = model.CSY;
    M.Q2Yhatcum = [];
    model = M;
    
    [~,nv]=size(model.W);
    switch prep
        case 'uv'
            WC = model.W';
        case 'pa'
            for k = 1:nv
                WC(k,:) = model.W(:,k)'./sqrt(model.StdX);
            end
        case 'mc'
            for k = 1:nv
                WC(k,:) = model.W(:,k)'./(model.StdX);
            end
    end
    Bool = exist('BT', 'var');
    if Bool == 1 && ~isempty(BT) && BT == 1 % bootstrapping
        
        f = @(x,y)o2pls_bootstrap(x,y,model.preprocessing,1,ncox,0,0);
        % % bootstrap
        [ci,bootstat] = bootci(500,{f,X,Y}, 'type', 'norm'); % bootstrap
        R = bootstat(:,1:size(X,2));
        median_R = median(R);
        [~,IX_R] = sort(abs(median_R),'descend');
        P = bootstat(:,size(X,2)+1:2*size(X,2)); % p-value
        median_P = median(P);
        W2_bt = bootstat(:,2*size(X,2)+1:3*size(X,2)); % W2
        median_W2_bt = median(W2_bt);
        if ~isempty(pCutoff)
            IX_P = find(median_P<pCutoff); %
        end
        
        
    else  % without bootstrapping
        [nsx,nvx]=size(X);
        buffer = repmat(model.T(:,1),[1,nvx]);
        CC = mean((buffer-repmat(mean(buffer),nsx,1)).*(X-repmat(mean(X),nsx,1)))/(nsx-1)*nsx;
        s1 = std(model.T(:,1));
        S1 = repmat(s1,[1,nvx]);
        S2 = std(X);
        median_R = CC./(S1.*S2);
        median_R(isnan(median_R))=0;
        median_P = tTest(median_R,nsx);
        if ~isempty(pCutoff)
            IX_P = find(median_P<pCutoff); %
        end
        
        switch prep
            case 'uv'
                %                 WC = model.W';
                W2=WC.*WC;
                for k = 1:nv
                    W2(k,:) = W2(k,:)/norm(W2(k,:), 'inf');
                end
            case 'pa'
                for k = 1:nv
                    %                     WC(k,:) = model.W(:,k)'./sqrt(model.StdX);
                    W2(k,:)=WC(k,:).*WC(k,:);
                    W2(k,:) = W2(k,:)/norm(W2(k,:), 'inf');
                end
            case 'mc'
                for k = 1:nv
                    %                     WC(k,:) = model.W(:,k)'./(model.StdX);
                    W2(k,:)=WC(k,:).*WC(k,:);
                    W2(k,:) = W2(k,:)/norm(W2(k,:), 'inf');
                end
        end
        median_W2_bt = W2(1,:);
    end
    
    Bool = exist('opt', 'var');
    if ~isempty(pCutoff) && ~isempty(r2cutoff)
        if Bool == 1 && ~isempty(opt) && opt == 0
            idx = find(median_W2_bt>r2cutoff);
        else
            II = find(median_W2_bt>r2cutoff);
            idx = intersect(IX_P, II);
        end
        
    else
        idx = [];
    end
    
    % % save output
    model.sig_idx = idx; % the significant ppm points
    model.R = median_R;
    model.P = median_P;
    %     model.WC_unscaled = WC.*model.StdX;
    model.W2 = median_W2_bt;
    model.WC = WC;
    model.signs = Signs;
    
    
else % CV
    
    % nrcv-fold cross validation
    
    block_num = floor(nsx/nrcv);
    %     Q2Yhatcum = zeros(1,nrcv);
    
    for cv = 1:nrcv % nrcv iterations of CV
        idx_test = nrcv*([1:block_num]-1)+cv;
        idx_tr = 1:nsx;
        idx_tr(idx_test) = [];
        X_test = X(idx_test,:);
        Y_test = Y(idx_test,:);
        X_tr = X(idx_tr,:);
        Y_tr = Y(idx_tr,:);
        
        [nsx_tr,nvx_tr]=size(X_tr);
        [nsy_tr,nvy_tr]=size(Y_tr);
        [nsx_test,nvx_test]=size(X_test);
        [nsy_test,nvy_test]=size(Y_test);
        switch prep
            
            case 'no'
                model.preprocessing='no';
                
            case 'mc'
                %                 disp('Meancentering')
                model.preprocessing='mc';
                X_tr=X_tr-repmat(mean(X_tr),nsx_tr,1);
                Y_tr=Y_tr-repmat(mean(Y_tr),nsy_tr,1);
                X_test=X_test-repmat(mean(X_test),nsx_test,1);
                Y_test=Y_test-repmat(mean(Y_test),nsy_test,1);
            case 'uv'
                % disp('Univariance scaling')
                model.preprocessing='uv';
                X_tr=X_tr-repmat(mean(X_tr),nsx_tr,1);
                Y_tr=Y_tr-repmat(mean(Y_tr),nsy_tr,1);
                X_tr=X_tr./repmat(std(X_tr),nsx_tr,1);
                Y_tr=Y_tr./repmat(std(Y_tr),nsy_tr,1);
                
                X_test=X_test-repmat(mean(X_test),nsx_test,1);
                Y_test=Y_test-repmat(mean(Y_test),nsy_test,1);
                X_test=X_test./repmat(std(X_test),nsx_test,1);
                Y_test=Y_test./repmat(std(Y_test),nsy_test,1);
            case 'pa'
                model.preprocessing='pa';
                X_tr=X_tr-repmat(mean(X_tr),nsx_tr,1);
                Y_tr=Y_tr-repmat(mean(Y_tr),nsy_tr,1);
                X_tr=X_tr./repmat(sqrt(std(X_tr)),nsx_tr,1);
                Y_tr=Y_tr./repmat(sqrt(std(Y_tr)),nsy_tr,1);
                
                X_test=X_test-repmat(mean(X_test),nsx_test,1);
                Y_test=Y_test-repmat(mean(Y_test),nsy_test,1);
                X_test=X_test./repmat(sqrt(std(X_test)),nsx_test,1);
                Y_test=Y_test./repmat(sqrt(std(Y_test)),nsy_test,1);
                
            otherwise
                disp('Unknown Preprocessing')
                model=[];
                return
                
        end
        
        % training
        model=O2pls(X_tr,Y_tr,nc,ncox,ncoy,'standard');
        if isempty(model)
            return;
        end
        %yhat
        modelPredy=O2plsPred(X_test, Y_test, model,ncox,ncoy,'x');
        %         %xhat
        %         modelPredx=O2plsPred(cvSet.xTest, cvSet.yTest, model,ioax-1,ioay-1,'y');
        
        % % the overall Q2
        SSY=sum(sum(Y_test.^2));
        Q2Yhatcum(cv) = 1-sum(sum((modelPredy.Yhat-Y_test).*(modelPredy.Yhat-Y_test)))/SSY;
    end
    model.Q2Yhatcum = mean(Q2Yhatcum(~isnan(Q2Yhatcum)));
    
    
    
end