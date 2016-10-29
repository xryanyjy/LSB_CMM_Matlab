function logProb = logNormalize(logProb, length)

max = -realmax('double');
for i = 1: length
    if logProb(i) > max
        max = logProb(i);
    end
end

sum = 0;
for i = 1:length
    logProb(i) = exp(logProb(i) - max);
    %disp(logProb);
    %disp(max);
    sum = sum + logProb(i);
end

%disp(sum);
%disp(logProb);

for i = 1:length
    logProb(i) = logProb(i) / sum;
   % disp(logProb);
end



