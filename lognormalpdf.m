function [h, hprime] = lognormalpdf(x, mu, sigma)

% Log of the univariate normal distribution with mean mu and
% standard deviation sigma.  Returns the log of the pdf and the
% derivative of the pdf at x.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt


h = 1/2*log(2*pi*sigma^2) -(x-mu).^2/(2*sigma^2);
hprime = -(x-mu)./sigma^2;