% stick breaking process with with logrithm values.  

%args:
% X: covariates of the data 
% W: weights for logistic function

%return value
% logTheta: each row is logrithm of multinomial probabilities of choosing each component

function logPhi = LSBCMM_logStickBreak(X,W)

link = X*W;
logV1 = -log(1+exp(-link));
for i = 1:size(logV1,1)*size(logV1,2)
    if -link(i) >30
        logV1(i) = link(i);
    end
end

logV0 = -log(1+exp(link));
for i = 1:size(logV0,1)*size(logV0,2)
    if link(i) >30
        logV0(i) = -link(i);
    end
end

accum = logV0;
accum(:,1) = 0;

for k = 2:size(logV1,2)
    accum(:,k) = accum(:,k-1) + logV0(:,k-1);
end
logPhi = logV1 + accum;
logPhi = LSBCMM_logNorm(logPhi);
