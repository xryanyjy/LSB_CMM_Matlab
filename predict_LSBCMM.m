function prob = predict_LSBCMM(object,testX)

testX_m = ones(size(testX,1),1);
testX = [testX_m, testX];
P = 1 ./ (1 + exp(-testX * object.W));
phi = LSBCMM_stickBreak(P);

prob = phi * (object.theta)';
