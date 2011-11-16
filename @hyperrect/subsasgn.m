function H = subsasgn(H, X, Q)
%SUBSASGN Indexed assignments for hyperrect objects
%
% H(I) = Q
%   assigns the values of Q into the elements of H specifed by the
%   subscript vector I.
%
% If H(I) = [], removes hyperrects at position I from the array (union).
%
% see also SUBSREF, HORZCAT
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)


if H.dims == 0
    return;
end

if numel(X) > 1,
    error('??? Attempt to reference field of non-structure array.');
else
    if (~strcmp(X.type,'()')),
        % only indexes in round brackets are allowed
        error(['Indexing with ''' X.type ''' not supported!']);
    end
end

idx = X.subs{1};

if isempty(Q) || Q.dims == 0
    H.L(:, idx) = [];
    H.H(:, idx) = [];
else
    assert(isa(Q, 'hyperrect'), 'The right hand side must be a hyperrect object.');
    assert(Q.dims == H.dims, 'Dimensions of both sides must agree.');
    
    H.L(:, idx) = Q.L;
    H.H(:, idx) = Q.H;
end

if isempty(H.L)
    H = hyperrect;
end

end
