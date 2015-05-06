function plotGmm(gmm, data, x)

% Plot a Gaussian mixture model described by the structure gmm.  This
% function works for 1D and 2D Gaussians.  GMM should have fields mu,
% s, and pi, where mu is a NxD matrix of the means of N Gaussians, s
% is a DxDxN matrix of precisions (inverse covariances), and pi is an
% N-vector of the relative weights of the Gaussians.  The second
% argument, DATA, is optional but should contain the data from which
% these Gaussians were generated, which will be plotted on the same
% graph.  The X argument, for 1D GMMs, is a 2-vector specifying the
% minimum and maximum X values over which to plot the Gaussians.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt

[N,D] = size(gmm.mu);

if(D == 1)
  % Number of standard deviations to go outside of the lowest and
  % highest means.
  xrange = 3;
  points = 222;
  
  sigmaSq = 1./gmm.s;
  sigma = sqrt(sigmaSq);
  
  if(nargin < 3)
    xmin = min(gmm.mu - xrange*sigma);
    xmax = max(gmm.mu + xrange*sigma);
    x = [xmin:(xmax-xmin)/points:xmax];
  end
  
  if(nargin > 1)
    N = numel(data);
    [histN, histX] = hist(data,round(N/10));
    dx = mean(diff(histX));
    bar(histX, histN/(N*dx), 1, 'w'), shading faceted
    hold on
  end
  
  for i=1:length(gmm.mu)
    y(i,:) = gmm.pi(i)/sqrt(2*pi*sigmaSq(i)) * ...
	exp(-(x-gmm.mu(i)).^2/(2*sigmaSq(i)));
  end
  tot = sum(y);

  plot(x, y, x, tot, '--');
  
  if(nargin > 1) hold off; end

elseif (D == 2)
  if(nargin > 1) 
      plot(data(:,1), data(:,2), '.k'); 
      hold on
  end
  w = gmm.pi / sum(gmm.pi);
  
  for i=1:size(gmm.mu,1);
    plotGauss2(gmm.mu(i,:), gmm.s(:,:,i), w(i));
    hold on
  end
  hold off


elseif (D == 3)
  grid on
  if(nargin > 1) 
      scatter3(data(:,1), data(:,2), data(:,3), 2, 'b'); 
      hold on
  end
  w = gmm.pi / sum(gmm.pi);
  
  for i=1:size(gmm.mu,1);
    plotGauss3(gmm.mu(i,:), gmm.s(:,:,i), w(i));
    hold on
  end
  hold off


else
  error('Can only plot 1D or 2D GMMs');
end
  


function plotGauss2(mu,S,w)
%plotGauss(mu,S)
%
%PLOT A 2D Gaussian
%This function plots the given 2D gaussian on the current plot.
 
t = -pi:.01:pi;
x = sin(t);
y = cos(t);
 
[vv,dd] = eig(S);
A = real((vv*sqrt(dd))');
z = [x' y']*A;

hand = plot(mu(1),mu(2),'X');
set(hand, 'Color', 1-[w w w], 'LineWidth', 3);
hand = plot(z(:,1)+mu(1),z(:,2)+mu(2));
set(hand, 'Color', 1-[w w w], 'LineWidth', 3);



function plotGauss3(mu, S, w)
% Plot a 3D Gaussian on the current plot

[u,s,v] = svd(S);
A = diag(1./diag(sqrt(s))) * v';
pts = 2*[A; -A] + repmat(mu, 6, 1);
for i=1:3
  h = plot3(pts(i+[0 3],1), pts(i+[0 3],2), pts(i+[0 3],3), 'r');
  set(h, 'linewidth', w*20);
end
