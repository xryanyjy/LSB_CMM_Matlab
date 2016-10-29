%fit a conditional mixture model
%input: 
% X: fitting covariates, with each row as an instance. 
% Y: target indicator matrix, with each row as 0-1 indicators of possible classes. One class is correct label.  
% sigma2: variance (\sigma^2) of the Gaussian prior over the weight matrix W. Play the role of C. 
% K: maximum number of components 
% alpha: prior of the Dirichlet process

%return value:
% model: the conditional mixture model

function model = lsbcmm_fit(X, Y, sigma2, K, alpha, marg_theta, quiet)

N = size(X,1);
temp_main = ones(N,1);
X = [temp_main, X];
D = size(X,2);
L = size(Y,2);

for i = 1:size(Y,1)
    k = 1;
    for j = 1:size(Y,2)
        if Y(i,j) == 1
            indY{i,1}(k) = j;
            k = k+1;
        end
    end
end

if ~quiet
    disp([num2str(L), ' classes']);
    disp(['Using ', num2str(K), ' components']);
end

model.likelihood = -Inf;
W = zeros(D,K);
W(1,:) = 0;
W(1,K) = Inf;

n_start = 1;
for trail = 1:n_start
    code = ovaCode(L,K);
    ealpha = code.ealpha;
    %in the case theta is not marginalized, theta is initialized as normalized ealpha instead of (ealpha - 1)
    colSums = sum(ealpha);
    for i = 1:size(ealpha,2)
        theta(:,i) = ealpha(:,i)./colSums(i);
    end
    
    logPhi = zeros(N,K);
    n_iter = 100;
    likelihood = zeros(n_iter,1);
    for iter = 1:n_iter
        %E step
        if marg_theta
            elog = ELogTheta(ealpha);
        else
            elog = log(theta);
            elog(isnan(elog)) = -10^9;
        end
        eVar = DecompE(logPhi',elog',indY,K,L);
        count = reshape(eVar.count,K,L)';
        ealpha = count + alpha;
        ephi = eVar.ephi';
        entropy = eVar.entropy;
        logpy = eVar.logpy;
        
        if marg_theta
            colSums1 = sum(ealpha);
            for col = 1:size(ealpha,2)
                for row = 1:size(ealpha,1)
                    theta(row,col) = ealpha(row,col)/colSums1(col);
                end
            end
        else
            %theta is the mode of the Dirichlet with ealpha. If alpha < 1, 
            %then the mode is not accurate.
            theta = ealpha - 1;
            for i = 1:length(theta)
                if theta(i) < 0
                    theta(i) = 1e-9;
                end
            end
            colSums2 = sum(ealpha);
            for col = 1:size(ealpha,2)
                for row = 1:size(ealpha,1)
                    theta(row,col) = ealpha(row,col)/colSums2(col);
                end
            end
        end
        %%
        % M step
        for k = 1:(K-1)   
             phik = ephi(:,k);
             phis = sum((ephi(:,(k+1):K)),2);
             Wk = W(:,k);
             Wk = lbfgs(@(Wk) ObjectivewithGrad(Wk,X,phis,phik,sigma2),Wk);
             W(:,k) = -Wk.X;
        end

        %%
        logPhi = LSBCMM_logStickBreak(X,W);
        % calculate the lower bound to decide when to stop
        if marg_theta
            matrix1 = ephi .* logPhi;
            term1 = sum(matrix1(:));
            term3 = logpy;
            term4 = entropy;
            elog = ELogTheta(ealpha);
            matrix2 = gammaln(sum(ealpha));
            sum_matrix2 = sum(matrix2(:));
            matrix3 = gammaln(ealpha);
            sum_matrix3 = sum(matrix3(:));
            matrix4 = (ealpha-alpha) .* elog;
            sum_matrix4 = sum(matrix4(:));
            term5 = sum_matrix2 - sum_matrix3 - K .* gammaln(alpha*L) + K .* (gammaln(alpha)*L) + sum_matrix4;
            matrix5 = W(2:D,1:(K-1)) .* W(2:D,1:(K-1));
            term6 = -0.5 * sum(matrix5(:)) ./ sigma2;
            likelihood(iter) = term1 + term3 + term4 - term5 +term6;
            if ~quiet
                disp(['Lower bound becomes ', num2str(likelihood(iter))]);
            end
        else
            % probability of W with Gaussian prior
            matrix5 = W(2:D,1:(K-1)).* W(2:D,1:(K-1));
            term1 = -0.5 * sum(matrix5(:)) / sigma2;
            
            % probability of \theta with Dirichlet prior
            matrix6 = log(theta);
            term2 = (alpha - 1) .* sum(matrix6(:));
            
            matrix7 = ephi .* logPhi;
            matrix8 = count .* log(theta);
            term3 = sum(matrix7(:)) + sum(matrix8(:));% Expectation of log p(y, z | X, W, \theta)
            term4 = entropy;% Entropy term 
            likelihood(iter) = term1 +term2 + term3 +term4;
            if ~quiet
                disp(['Lower bound becomes', num2str(likelihood(iter))]);
            end
        end
        if iter > 10
            matrix9 = likelihood((iter-2):iter);
            curLike = mean(matrix9(:));
            matrix10 = likelihood((iter-5):(iter-3));
            lastLike = mean(matrix10(:));
            if (curLike-lastLike)<1e-4 * abs(curLike)
                break;
            end
        end
    end
    if model.likelihood < likelihood(n_iter)
        model.W = W;
        model.theta = theta;
        model.likelihood = likelihood(n_iter);
    end
end






