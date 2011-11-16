function obj = removethin(obj)
%REMOVETHIN Remove thin sets from a hyper-rectangle array.
%      B = removethin(A)
%
%   A is an array of hyper-rectangles (or single one).
%   All thin sets (i.e. l_i = h_i for some i) are removed from A. The
%   result is saved in B.  If all sets in A are thin, then B is the empty
%   hyper-rectangle.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if obj.dims == 0, return; end

thinIdx = any(obj.L >= obj.H);  % Indices of thin sets

if all(thinIdx)
    % All are thin
    obj = hyperrect;
else
    % Only keep those that are full-dim
    obj.L(:, thinIdx) = [];
    obj.H(:, thinIdx) = [];
end

end

