function iss = ge(A, B)
%GE Check if a hyperrectangle set is a superset of another.
%
%   H1 >= H2
%       Returns true if H2 is a subset of H1; both of them can be an array
%       (union of hyper-rectangles).
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

iss = le(B, A);

end
