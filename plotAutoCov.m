function plotAutoCov(Samp)

% Plot the autocovariances of the various parameters as a function
% of lag in the sampling.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt


acK = ifftshift(xcov([Samp.k], 'coeff'));
% $$$ acR = ifftshift(xcov([Samp.r], 'coeff'));
% $$$ acW = ifftshift(xcov([Samp.w], 'coeff'));
acLambda = ifftshift(xcov([Samp.lambda], 'coeff'));
acAlpha = ifftshift(xcov([Samp.alpha], 'coeff'));
acBeta = ifftshift(xcov([Samp.beta], 'coeff'));

x = [1:1000];

plot(x, acK(x), x, acAlpha(x), x, acLambda(x), x, ...
    acBeta(x), [1 1000], [0 0]);
% $$$ plot(x, acK(x), x, acAlpha(x), x, acW(x), x, acLambda(x), x, ...
% $$$     acR(x), x, acBeta(x), [1 1000], [0 0]);

legend('k', '\alpha', '\lambda', '\beta');
% $$$ legend('k', '\alpha', 'w', '\lambda', 'r', '\beta');
title('Autocovariances of Hyperparameters');
xlabel('Lag')