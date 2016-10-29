% stick breaking process with probability P

% input:
% P: each row are Bernoulli probabilities of choosing each component

% return value:
% phi: each row are multinomial probabilities of choosing each component
function phi = LSBCMM_stickBreak(P)

phi = P;
temp = 1-P;

for k = 2:size(P,2)
    for i = 1:size(phi,1)
        phi(i,k) = temp(i,k-1) * P(i,k);
        temp(i,k) = temp(i,k-1) * (1-P(i,k));
    end
end

 