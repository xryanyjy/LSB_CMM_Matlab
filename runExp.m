% An example using this algorithm to run the Lost dataset
load lost_reduce
target = target';
partial_target = partial_target';
data = zscore(data);
acc = zeros(10,1);
 
 for i = 1:10
     train_data = data(tr_idx{i,1},:);
     train_target = partial_target(tr_idx{i,1},:);
     test_data = data(te_idx{i,1},:);
     test_target = target(te_idx{i,1},:);
 
     model = lsbcmm_fit(train_data,train_target,1,80,0.05,true,false);
     prob = predict_LSBCMM(model,test_data);

     acc(i) = accuracy(prob,test_target);
     acc_mean = mean(acc);
      acc_std = std(acc);
end
 


