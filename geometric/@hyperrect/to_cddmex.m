function S = to_cddmex(H)
%TO_CDDMEX Convert to the H-form structure of cddmex.
%   S = to_cddmex(H)
%
%See also: hyperrect/to_polytope
%
% (C) 2012 by Truong X. Nghiem (nghiem@seas.upenn.edu)

[Ai, Bi, Ae, Be] = getconstraints(H);
S = struct('A', [Ae; Ai], 'B', [Be; Bi], 'lin', 1:length(Be));

end

