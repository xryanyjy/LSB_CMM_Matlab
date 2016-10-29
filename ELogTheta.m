% calculate expectation of log theta and theta is from dirichlet distribution with alpha as parameter
% input:
% ealpha: parameter of the Dirichlet distribution

% return value:
% mean.log: expected value of log(\theta_i) for each i

function mean_log = ELogTheta(ealpha)
mean_log = ealpha;
for k = 1:size(ealpha,2)
  %disp(ealpha(:,k)');
  x = drchrnd(ealpha(:,k)',10000);
  %disp(x);
  mean_log(:,k) = mean(log(x))';
  %disp(mean_log);
  if any(mean_log == -Inf)
      error('The expectation becomes minus infinity');
  end
end