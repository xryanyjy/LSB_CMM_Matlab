function acc = accuracy(prob,test_target)
[~,index_model] = max(prob,[],2);
%[max_b,index_num] = max(test_target,[],2);
sum = 0;
test_target = full(test_target);
for i = 1:size(test_target,1)
%     if index_model(i) == test_target(i);
    if test_target(i,index_model(i))==1 
        sum = sum + 1;
    end
end
acc = sum / size(test_target,1);
