% the objective value of M step
function [f,g] = ObjectivewithGrad(W,X,phik,phis,sigma2)

f = (W(2:end)'*W(2:end)/2)/sigma2;
value = X*W;
f = f + (phik+phis)'*log(1+exp(value)) - phik'*value;
if nargout > 1
g=W / sigma2;
g(1) =0; 
g = g+ X'*(exp(value)./(1+exp(value)).*(phik+phis)-phik);
end

end







