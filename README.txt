Michael Mandel
CS 4771 Final Project
The Infinite Gaussian Mixture Model
Prof. Tony Jebara
May 5, 2005

In order to generate the test data used in the paper, just make this
call in matlab:
[Y,z] = drawGmm([-3 3], [1 10], [1 2], 500);

In order to run the infinite GMM on the data for 10000 iterations,
make this call:
Samp = igmm_uv(Y, 10000);

It's as easy as that.

If you want to run the regular univariate Gibbs Sampler on the data,
do this: 
[mu,sigSq,p,z,churn] = gibbsGmm(Y,2,0,100,2,1,2,1000);

The igmm for multivariate data is in igmm_mv.m, which uses
logmvbetpdf.m instead of the logbetapdf.m used by igmm_uv.m.
Otherwise, both igmms are self-contained.

To generate multivariate data, use e.g.
S = [2 1; 1 2]; S(:,:,2) = S
[Y,z] = drawGmm([3 -3; -3 3], S, [1 1], 100);
Samp = igmm_mv(Y, 10000);

To generate figure 3 in the paper, use the function plotAutoCov.m



=====================================================
COPYRIGHT / LICENSE
=====================================================
All code was written by Michael Mandel, and is copyrighted under the
(lesser) GPL:
  Copyright (C) 2005  Michael Mandel
 
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; version 2.1 or later.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 
The authors may be contacted via email at: mim at ee columbia edu
