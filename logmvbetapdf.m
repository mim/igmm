function [h, hprime] = logmvbetapdf(y, S, W)

% Generate the log of the pdf for y=beta-D+1 in rasmussen's infinite
% gmm, eq 9 for the multivariate case.  The i-th covariance matrix is
% s(:,:,i).  To get samples from beta, take samples from this
% function and add D-1 to them.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt


if(y < 0) y, end
k = size(S,3);
d = size(S,1);
N = length(y);

ldW = log(det(W));
for i=1:k
  ldS(i) = log(det(S(:,:,i)));
  trWS(i) = trace(W*S(:,:,i));
end
yd1 = y+d-1;

sumGammaLn = 0; sumPsi = 0;
for i=0:d-1
  sumGammaLn = sumGammaLn + gammaln(.5*(y+i));
  sumPsi = sumPsi + psi(.5*(y+i));
end

if(k < N)
  sum_ldS_trWS = zeros(size(y));
  for i=1:k
    sum_ldS_trWS = sum_ldS_trWS + (y+2)/2*ldS(i) - yd1/2*trWS(i);
  end
else
  for i=1:N
    sum_ldS_trWS(i) = sum((y(i)+2)/2*ldS - yd1(i)/2*trWS);
  end
end

h = -3/2*log(yd1) - d./(2*yd1) - k*(y*d*log(2)/2 + sumGammaLn) + ...
    yd1*k/2.*(d*log(yd1) + ldW) + sum_ldS_trWS;
hprime = -3./(2*yd1) + d./(2*yd1.^2) - k*(d*log(2)/2 + sumPsi/2) + ...
    k/2*(d*log(yd1) + ldW + d) + .5*sum(ldS - trWS);
