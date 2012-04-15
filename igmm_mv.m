function Samp = igmm_mv(Y, Nsamp, init, debg)

% Samp = igmm_mv(Y, Nsamp, init)
%
% Implementation of Rasmussen's infinite GMM for multivariate data.
% Use Gibbs sampling to draw Nsamp samples from the posterior
% distribution of an infinite GMM given some data Y, where each row is
% a data point.  Based on a heirarchical graphical model where
% hyperparameters are also inferred, so there are no parameters to
% tweak by hand.  Samp is a vector of structures, where each element
% has fields mu, s, lambda, r, beta, w, alpha, and k as described in
% the paper.  Start the Gibbs sampler with the sample INIT, if
% supplied.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt


if(nargin < 4) debg = 0; end

[N,D] = size(Y);
mu_y = mean(Y);
sigSq_y = cov(Y);
sigSqi_y = inv(sigSq_y);
arsStart = [.2 2];

if(nargin < 3 || isempty(init))
  % Start off with one class
  c = ones(1,N);
  Samp(1).k = 1;
  Samp(1).mu = mu_y;
  Samp(1).s = sigSqi_y;
  Samp(1).lambda = drawNormal(mu_y, sigSq_y);
  Samp(1).r = drawWishart(D, sigSqi_y);
  Samp(1).beta = 1/drawGamma(1,1/D) + D-1;
  Samp(1).w = drawWishart(D, sigSq_y);
  Samp(1).alpha = 1/drawGamma(1,1);
  Samp(1).pi = 1;
  Samp(1).Ic = 1;
else
  % Initialize with supplied sample
  Samp(1) = init;
  
  % Since c is not included in the intialization, sample it
  c = drawMultinom(normalLikeSame(Y, init.mu, init.s));
end

Ic = 1;

% Go!
for i=2:Nsamp
  prt(debg, 1, '########### ', i);
  
  % Make aliases for more readable code
  k = Samp(i-1).k;    mu = Samp(i-1).mu;    s = Samp(i-1).s;
  beta = Samp(i-1).beta;    r = Samp(i-1).r;    w = Samp(i-1).w;
  lambda = Samp(i-1).lambda;    alpha = Samp(i-1).alpha;
  
  % Find the popularity of the classes
  nij = tabulate(c); nij = nij(:,2);
  prt(debg, 2, 'nij = ', nij')
  
  %%%%% Mu
  for j=1:k
    inClass = find(c == j);
    n = numel(inClass);
    if(n <= 0) ybar = 0;  
    else       ybar = mean(Y(inClass,:),1);
    end

    tmp_sigSq = inv(n*s(:,:,j) + r);
    tmp_mu = (n*ybar*s(:,:,j) + lambda*r)*tmp_sigSq;
    Samp(i).mu(j,:) = drawNormal(tmp_mu, tmp_sigSq);
  end
  prt(debg, 3, 'mu = ', Samp(i).mu);

  
  %%%%% Lambda
  tmp_sigSq = inv(sigSqi_y + k*r);
  tmp_mu = (mu_y*sigSqi_y + sum(mu,1)*r) * tmp_sigSq;
  Samp(i).lambda = drawNormal(tmp_mu, tmp_sigSq);
  prt(debg, 3, 'lambda = ', Samp(i).lambda);
  

  %%%%%% R
  mMinL = mu - repmat(lambda, k, 1);
  Samp(i).r = drawWishart(k+1, 1/(k+1)*inv(sigSq_y + mMinL'*mMinL));
  prt(debg, 3, 'r = ', Samp(i).r);

  
  %%%%%% S
  for j=1:k
    inClass = find(c == j);
    n = numel(inClass);
    if(n <= 0) sbar = 0;
    else       
      yMinMu = Y(inClass,:) - repmat(mu(j,:), n, 1);
      sbar = yMinMu' * yMinMu;
    end

    tmp_a = beta+nij(j);
    tmp_b = tmp_a*inv(w*beta + sbar);
    Samp(i).s(:,:,j) = drawWishart(tmp_a, tmp_b);
  end
  prt(debg, 3, 's = ', Samp(i).s);

  
  %%%%%% W
  tmp_a = k*beta+1;
  tmp_b = tmp_a * inv(sigSqi_y + beta*sum(s,3));
  Samp(i).w = drawWishart(tmp_a, tmp_b);
  prt(debg, 3, 'w = ', Samp(i).w);
  
  
  %%%%%% Beta
  Samp(i).beta = ars(@logmvbetapdf, {s,w}, 1, arsStart, [0 inf])+D-1;
  prt(debg, 2, 'beta = ', Samp(i).beta)
  
  
  %%%%%% Alpha
  Samp(i).alpha = ars(@logalphapdf, {k, N}, 1, arsStart, [0 inf]);
  prt(debg, 2, 'alpha = ', Samp(i).alpha);
  
  
  %%%%%% C
  % Samples from priors, which could be swapped in if we need a new
  % class.  Only needed for i>=Ic
  mu_prop = [zeros(Ic-1,D); drawNormal(lambda, pinv(r), N-Ic+1)];
  s_prop = ones(D,D,N);
  wi = inv(w);
  [s_prop(:,:,Ic), WiCh] = drawWishart(beta, wi);
  for prop = Ic+1:N
    s_prop(:,:,prop) = drawWishart(beta, wi, WiCh);
  end

  % find the liklihoods under samples from the prior for i>=Ic
  unrep_like = alpha/(N-1+alpha) * normalLikeDiff(Y, mu_prop, s_prop);

  % Find the likelihoods of the observations under the *new* gaussians
  mu = Samp(i).mu;   s = Samp(i).s;
  rep_like = normalLikeSame(Y, mu, s);
  
  % Calculate the priors, specific to each datapoint, because
  % counting over everyone else
  pri = repmat(nij/(N-1+alpha), 1, N);
  idxs = sub2ind(size(pri), c, [1:N]);
  pri(idxs) = pri(idxs) - 1/(N-1+alpha);
  
  % Tweak probabilities for classes with one member
  for cl=find(nij == 1)'
    idx = find(c == cl);
    unrep_like(idx) = 0;
    pri(cl,idx) = pri(cl,idx) + 1/(N-1+alpha);
  end
  
  % Assign datapoints to classes, get rid of empty classes, add new
  % ones
  like = [pri .* rep_like; unrep_like];
  prt(debg, 2, 'like = ', like);
  cn = drawMultinom(like);
  [c,keep,to_add,Ic] = renumber(cn, c, k, Ic);
  Samp(i).mu = [mu(keep,:); mu_prop(to_add,:)];
  Samp(i).s = s(:,:,keep);
  Samp(i).s(:,:,end+1:end+numel(to_add)) = s_prop(:,:,to_add);
  Samp(i).k = size(Samp(i).mu, 1);

  prt(debg, 4, 'mu = ', Samp(i).mu);
  prt(debg, 4, 'k = ', Samp(i).k);

  
  % Find relative weights of components
  nij = tabulate(c); nij = nij(:,2);
  Samp(i).pi = nij / sum(nij);

  Samp(i).Ic = Ic;
  prt(debg, 2, 'Ic = ', Ic);
end



function x = drawNormal(mu, sigSq, N)
% Draw N samples from a Gaussian with mean mu and covariance sigSq.
% Mu is a row vector.
if(nargin < 3) N = 1; end
D = size(mu,2);
[u,s,v] = svd(sigSq);
sig = sqrt(s)*v';
x = randn(N,D)*sig + repmat(mu, N, 1);


function pr = normalLikeSame(y, mu, s)
% Evaluate the likelihood of data y under each of the gaussians
% with mean mu(j,:) and precision s(:,:,j).
[Ny, D] = size(y);
Nmu = size(mu,1);
for j=1:Nmu
  yMinMu = (y - repmat(mu(j,:), Ny, 1))';
  pr(j,:) = sqrt(det(s(:,:,j)/(2*pi))) * ...
      exp( -1/2 * sum( yMinMu .* (s(:,:,j) * yMinMu) ) );
end


function pr = normalLikeDiff(y, mu, s)
% Evaluate the likelihood of data point y(j,:) under the gaussian
% with mean mu(j,:) and precision s(:,:,j).
[Nmu, D] = size(mu);
for j=1:Nmu
  yMinMu = y(j,:) - mu(j,:);
  pr(j) = sqrt(det(s(:,:,j)/(2*pi))) * exp(-1/2*yMinMu * s(:,:,j) * yMinMu');
end


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
x = gamrnd(shape, mean./shape);



function [X,Mch] = drawWishart(shape, mean, varargin)
% Draw a wishart-distributed random matrix having shape and mean
% parameters given by the arguments.
[X,Mch] = wishrnd(mean/shape, shape, varargin{:});
while(det(X) <= 0) 
% $$$   warning('Wishart drew nonpositive definite matrix'); 
  X = wishrnd(mean/shape, shape, Mch);
end



function prt(debg, level, txt, num)
% Print text and number to screen if debug is enabled.
if(debg >= level)
  if(numel(num) == 1)
    disp([txt num2str(num)]);
  else
    disp(txt)
    disp(num)
  end
end
