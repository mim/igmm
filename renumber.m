function [c,keep,to_add,Ic] = renumber(cn, c, k, Ic)
% Assign new classes to points.  From Ic up to the first addition
% of a new class, c = cn, after the addition of the new class, c =
% c, and Ic is set to the index of the first point that wasn't
% relabeled.  A new class is indicated by cn = k+1.  If there is
% no new class added between Ic and the end of cn, reset Ic to 1.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt


Icn = find(cn(Ic:end) == k+1)+Ic-1;
is_new = 1-isempty(Icn);

if(is_new)
  Icn = Icn(1);
  c(Ic:Icn) = cn(Ic:Icn);
  Ic = Icn+1;
  to_add = Icn;
else
  c(Ic:end) = cn(Ic:end);
  to_add = [];
  Ic = 1;
end

% count up how many members of each class
tab = tabulate(c);

% get rid of empty classes
keep = find(tab(:,2) ~= 0);
tab = tab(keep,:);

% don't include the extra class if there is one
keep = keep(1:end-is_new);

% count the number of classes to keep
kn = size(tab,1);

% relabel points so that each group has at least one point in it
for j=1:kn
  c(find(c == tab(j,1))) = j;
end
