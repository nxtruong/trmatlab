function f = isfulldim(obj)
%ISFULLDIM Check if a hyper-rectangle is fully dimensional.
%      f = isfulldim(H)
%           Returns true if H is fully dimensional.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

f = (obj.dims > 0) && any(all(obj.L < obj.H));

end