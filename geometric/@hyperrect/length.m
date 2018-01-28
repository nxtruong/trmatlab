function n = length(obj)
%LENGTH   Number of hyperrectangles in the union.
%   LENGTH(H) returns the number of hyperrectangles in the union H.
%   Union of hyperrectangles can be constructed using the standard
%   concatenation operations, e.g.
%       H = [H1, H2, H3]
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if obj.dims == 0
    n = 0;
else
    n = size(obj.L, 2);
end

end

