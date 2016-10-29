%Normalize each row (logrithm values) to be logrithm of a multinomial probability. 
%input
%logProb: each row is logrithm of values proportional to some probability 
%return value
%logProb: each row is logrithm of multinomial probability 
function logProb = LSBCMM_logNorm(logProb)
%------------------------------------------
mv = max(logProb,[],2);
%t_logProb = logProb';
%mean_value = mean(logProb);
logProb_col = size(logProb,2);%列数
logProb_row = size(logProb,1);%行数
for i=1:logProb_row
    for j = 1:logProb_col
        logProb(i,j) = logProb(i,j) - mv(i);
    end
end
%logProb = t_logProb';
prob = exp(logProb);
norm = sum(prob,2);
%lognorm = matrix(log(norm),ncol=1);
%lognorm = zeros(length(norm),1);
lognorm = log(norm);

for i = 1:size(lognorm,1)
    for j = 1:size(logProb,2)
        rep_logProb(i,j) = lognorm(i);
    end
end
logProb = logProb - rep_logProb;









