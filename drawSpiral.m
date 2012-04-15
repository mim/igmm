function Y = drawSpiral(N, std)

% Y = drawSpiral(N, std)
%
% Draw N points from a noisy 3D spiral similar to the one used in
% Rasmussen's paper, which was taken from Ueda et al (1998).  Std is
% the standard deviation of noise around the spiral.  The default
% parameters should give something like Ueda's spiral.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt


if(nargin < 2) std = .05; end
if(nargin < 1) N = 800; end

t = rand(1,N)*4*pi + 2*pi;
Y = [10*cos(t)./t;-10*sin(t)./t; t/(4*pi)]';
Y = Y + randn(size(Y))*std;
