function [mu, sigmaSq, p, z, churn] = ...
    gibbsGmm(Y, k, m, etaSq, nu0, nu0lambda0, alpha, Nsamp)

% Use Markov chain Monte Carlo simulation to cluster the data Y into a
% mixture of k univariate Gaussians.  Priors on variables are: mu ~
% N(m, etaSq), sigmaSq ~ Wishart(nu0, lambda0), pi ~ dirichlet(alpha/k).
% Outputs of function are samples from the posterior distributions, so
% that theta(i) = [mu(i,:) sigma(i,:) z ], i = 1..Nsamp

N = length(Y);

% Randomly assign to classes, initialize stats to the means and
% vars of those classes.
z = drawMultinom(ones(k,N));
for j=1:k
  yj = Y(find(z == j));
  mu(1,j) = yj(unidrnd(numel(yj)));
  sigmaSq(1,j) = std(yj).^2;
end
p(1,:) = full(sparse(1, z, 1, 1, k));

% Go!
for i=2:Nsamp
  % Mu
  for j=1:k
    n = sum(z == j);
    if(n <= 0) ybar = 0;  
    else       ybar = mean(Y(find(z == j)));
    end

    tmp_sigSq = 1/(n/sigmaSq(i-1,j) + 1/etaSq);
    tmp_mu = tmp_sigSq*(n*ybar/sigmaSq(i-1,j) + m/etaSq);
    mu(i,j) = drawNormal(tmp_mu, tmp_sigSq);
  end

  % Sigma
  for j=1:k
    inClass = z == j;
    n = sum(inClass);
    if(n <= 0) sigbar = 0;
    else       sigbar = sum((Y(find(inClass)) - mu(i-1,j)).^2);
    end
    
    tmp_nu = nu0+n;
    tmp_nu_lambda = (nu0lambda0 + sigbar);
    sigmaSq(i,j) = drawInvChiSq(tmp_nu, tmp_nu_lambda);
  end
  
  % z \in {1..k}
  for j=1:k
    tmp_pr(j,:) = normalLike(Y, mu(i-1,j), sigmaSq(i-1,j));
  end
  n = tabulate(z);
  n = n(:,2)';
% $$$   n = full(sparse(1, z, 1, 1, k));
  
  % Scale likelihoods by class memberships times prior
  pri = repmat((n'+alpha/k)/(sum(n)-1+alpha), 1, N);
  idxs = sub2ind(size(pri), z, [1:N]);
  pri(idxs) = pri(idxs) - 1/(sum(n)-1+alpha);
  
  tz = drawMultinom(pri .* tmp_pr);
  churn(i) = sum(tz ~= z);
  z = tz;
  p(i,:) = n;
  
% $$$   plotGmm(mu(i,1), mu(i,2), sigmaSq(i,1), sigmaSq(i,2), p(i));
% $$$   pause(.1)
end


function x = drawNormal(mu, sigSq)
% Draw one sample from a Gaussian with mean mu and variance sigSq
x = randn(1)*sqrt(sigSq) + mu;


function pr = normalLike(y, mu, sigSq)
% Evaluate the likelihood of the points y under the Gaussian with mean
% mu and variance sigSq
pr = 1/sqrt(2*pi*sigSq) .* exp(-(y-mu).^2/(2*sigSq));


function x = drawInvChiSq(nu, nu_lambda)
% Draw one sample from an inverse chi square distribution with
% parameters nu and lambda
x = nu_lambda / chi2rnd(nu);


function x = drawBeta(a, b)
% Draw one sample from a Beta distribution with parameters a and b
x = betarnd(a,b);


function x = drawBernoulli(p)
% Draw bernoulli random variables with probability p of getting 1.
% x is the same size as p.
x = rand(size(p)) < p;
