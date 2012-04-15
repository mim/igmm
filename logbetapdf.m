function [h, hprime] = logbetapdf(beta, s, w)

% Generate the log of the pdf for the beta variable in rasmussen's
% infinite gmm, eq 9.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt


if(beta < 0)
  beta
end
k = length(s);
h = -k * gammaln(beta/2) - 1./(2*beta) + (k*beta-3)./2.*log(beta/2) + ...
    (beta/2)*sum(log(w) + log(s)) - beta*sum(s*w/2);
hprime = -k/2*psi(beta/2) + 1./(2*beta.^2) + k/2.*log(beta/2) + ...
    (k*beta-3)./(2*beta) + sum(log(s) + log(w) - s*w)/2;
