function Samp = igmm_uv(Y, Nsamp, debg)

% Implementation of Rasmussen's infinite GMM for univariate data.  Use
% Gibbs sampling to draw Nsamp samples from the posterior distribution
% of an infinite GMM given some data Y.  Based on a heirarchical
% graphical model where hyperparameters are also inferred, so there
% are no parameters to tweak by hand.  Samp is a vector of structures,
% where each element has fields mu, s, lambda, r, beta, w, alpha, and
% k as described in the paper.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt

% TODO:
% Handle case where n_{-i,j}=0, but c_i=j

if(nargin < 3) debg = 0; end

N = length(Y);
mu_y = mean(Y);
sigSq_y = var(Y);
sigSqi_y = 1/sigSq_y;
arsStart = [.1 4];

% Start off with one class
c = ones(1,N);
Samp(1).k = 1;
Samp(1).mu = mu_y;
Samp(1).s = sigSqi_y;
Samp(1).lambda = drawNormal(mu_y, sigSq_y);
Samp(1).r = drawGamma(1, sigSqi_y);
Samp(1).beta = 1/gamrnd(1,1);
Samp(1).w = drawGamma(1, sigSq_y);
Samp(1).alpha = 1/gamrnd(1,1);
Samp(1).pi = 1;
Samp(1).Ic = 1;

Ic = 1;

% Go!
for i=2:Nsamp
  prt(debg, '########### ', i);
  
  % Make aliases for more readable code
  k = Samp(i-1).k;    mu = Samp(i-1).mu;    s = Samp(i-1).s;
  beta = Samp(i-1).beta;    r = Samp(i-1).r;    w = Samp(i-1).w;
  lambda = Samp(i-1).lambda;    alpha = Samp(i-1).alpha;
  
  % Find the popularity of the classes
  nij = tabulate(c); nij = nij(:,2);
  prt(debg, 'nij = ', nij')
  
  %%%%% Mu
  for j=1:k
    inClass = find(c == j);
    n = numel(inClass);
    if(n <= 0) ybar = 0;  
    else       ybar = mean(Y(inClass));
    end

    tmp_sigSq = 1/(n*s(j) + r);
    tmp_mu = tmp_sigSq*(n*ybar*s(j) + lambda*r);
    Samp(i).mu(j) = drawNormal(tmp_mu, tmp_sigSq);
  end
  prt(debg, 'mu = ', Samp(i).mu);

  
  %%%%% Lambda
  tmp_sigSq = 1/(sigSqi_y + k*r);
  tmp_mu = tmp_sigSq*(mu_y*sigSqi_y + r*sum(mu));
  Samp(i).lambda = drawNormal(tmp_mu, tmp_sigSq);
  prt(debg, 'lambda = ', Samp(i).lambda);
  

  %%%%%% R
  Samp(i).r = drawGamma(k+1, (k+1)/(sigSq_y + sum((mu-lambda).^2)));
  prt(debg, 'r = ', Samp(i).r);

  
  %%%%%% S
  tmp_a = []; tmp_b = [];
  for j=1:k
    inClass = find(c == j);
    n = numel(inClass);
    if(n <= 0) sbar = 0;
    else       sbar = sum((Y(inClass) - mu(j)).^2);
    end

    tmp_a(j) = beta+nij(j);
    tmp_b(j) = tmp_a(j) / (w*beta + sbar);
  end
  Samp(i).s = drawGamma(tmp_a, tmp_b);
  prt(debg, 's = ', Samp(i).s);

  
  %%%%%% W
  Samp(i).w = drawGamma(k*beta+1, (k*beta+1)/(sigSqi_y + beta*sum(s)));
  prt(debg, 'w = ', Samp(i).w);
  
  
  %%%%%% Beta
  % is beta log concave or log(beta)?  Looks like beta itself is.
  Samp(i).beta = ars(@logbetapdf, {s, w}, 1, arsStart, [0 inf]);
  prt(debg, 'beta = ', Samp(i).beta)
  
  
  %%%%%% Alpha
  % is alpha log concave or log(alpha)? looks like alpha itself is.
  Samp(i).alpha = ars(@logalphapdf, {k, N}, 1, arsStart, [0 inf]);
  prt(debg, 'alpha = ', Samp(i).alpha);
  
  
  %%%%%% C
  % Samples from priors, which could be swapped in if we need a new
  % class.  Only needed for i>=Ic
  mu_prop = [zeros(1,Ic-1) drawNormal(lambda, 1/r, N-Ic+1)];
  s_prop = [ones(1,Ic-1) drawGamma(beta*ones(1,N-Ic+1), 1/w)];

  % Find the likelihoods of the observations under the *new*
  % gaussians for i>=Ic
  unrep_like = alpha/(N-1+alpha) * normalLike(Y, mu_prop, s_prop);

  mu = Samp(i).mu;   s = Samp(i).s;
  rep_like = [];
  for j=1:k
    rep_like(j,:) = normalLike(Y, mu(j), s(j));
  end
  
  % Calculate the priors, specific to each datapoint, because
  % counting over everyone else
  pri = repmat(nij/(N-1+alpha), 1, N);
  idxs = sub2ind(size(pri), c, [1:N]);
  pri(idxs) = pri(idxs) - 1/(N-1+alpha);
  
  % Assign datapoints to classes, get rid of empty classes, add new
  % ones
  cn = drawMultinom([pri .* rep_like; unrep_like]);
  [c,keep,to_add,Ic] = renumber(cn, c, k, Ic);
  Samp(i).mu = [mu(keep) mu_prop(to_add)];
  Samp(i).s = [s(keep) s_prop(to_add)];
  Samp(i).k = length(Samp(i).mu);

  prt(debg, 'mu = ', Samp(i).mu);
  prt(debg, 'k = ', Samp(i).k);

  
  % Find relative weights of components
  nij = tabulate(c); nij = nij(:,2);
  Samp(i).pi = nij / sum(nij);

  Samp(i).Ic = Ic;
end



function x = drawNormal(mu, sigSq, N)
% Draw one sample from a Gaussian with mean mu and variance sigSq
if(nargin < 3) N = 1; end
x = randn(1,N)*sqrt(sigSq) + mu;


function pr = normalLike(y, mu, s)
% The likelihood of data y under the gaussian with mean mu and inverse
% variance s.  If y, mu, and s are the same size, evaluate the
% likelihood of each point in y under a different mu and s, if mu
% and s are scalars, evaluate all points under the same gaussian.
pr = sqrt(s/(2*pi)) .* exp(-s.*(y-mu).^2/2);


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


function x = drawGamma(shape, mean)
% Draw a gamma-distributed random variable having shape and mean
% parameters given by the arguments.  Translate's Rasmussen's shape
% and mean notation to mathworld's and mathworks' alpha and theta
% notation.  When rasmussen writes G(beta, w^-1), matlab expects
% G(beta, w^-1/beta).
x = gamrnd(shape/2, 2*mean./shape);


function prt(debg, txt, num)
% Print text and number to screen if debug is enabled.
if(debg)
  disp([msg num2str(num)]);
end
