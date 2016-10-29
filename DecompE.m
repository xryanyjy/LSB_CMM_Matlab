function mAns = DecompE(logPhi, elog, Y, K, L)

N = numel(Y);
mAns.logpy = 0;
mAns.entropy = 0;
mAns.ephi = zeros(K,N);
mAns.count = zeros(K*L,1);

% for i = 1:L
%     fullLabels(i) = i;
% end
insEphi = zeros(K*L,1);

for n = 1:N
        rY = Y{n,1};
        %disp(n);
        insLab = rY;
        nLab = numel(rY);
   
    
    for c = 1:nLab
        lab = insLab(c) - 1;
        %if lab == -1
        %    disp(insLab(c));
        %    disp(c);
        %end
        %disp(lab);
        for k = 1:K
            insEphi((c-1) * K +k) = logPhi((n-1) * K + k) + elog(lab * K +k);
            %disp(n);
        end
    end
    %disp(insEphi);
    %disp(nLab*K);
    insEphi = logNormalize(insEphi,nLab * K);
    %disp(insEphi);
    %ephi has been initialized at the beginning
    for c = 1:nLab
        lab = insLab(c) - 1;
        for k = 1:K
            %disp(insEphi);
            mAns.ephi((n-1) * K +k) = mAns.ephi((n-1) * K + k) + insEphi((c-1) * K + k);
            mAns.count(lab * K + k) = mAns.count(lab * K + k) + insEphi((c-1) * K + k);
            %disp(mAns.count);
        
            if insEphi((c-1) * K + k) > 0
                mAns.entropy = mAns.entropy - insEphi((c-1) * K + k) * log(insEphi((c-1) * K + k));
                mAns.logpy = mAns.logpy + insEphi((c-1) * K + k) * elog(lab * K + k);
            end
        end
    end
end
clear insEphi;

    
    
    
    
    
    
    
    