function b = isbounded(obj)
%ISBOUNDED Check if a (union of) hyper-rectangle is bounded.
%      b = isbounded(H)
%           Returns true if H is bounded.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

b = (obj.dims == 0) || ...
    (~any(any(isinf(obj.L))) && ~any(any(isinf(obj.H))));

end

