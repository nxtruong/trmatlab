function H = subsref(H, X)
%SUBSREF Indexed referencing for hyperrect objects.
%
% H(I) is an array formed from the elements of H specifed by the
%      subscript vector I.
%
% Indexing by logicals is supported as well.
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

if numel(X) > 1
    error('??? Attempt to reference field of non-structure array.');
else
    if (~strcmp(X.type,'()'))
        % only indexes in round brackets are allowed
        error(['Indexing with ''' X.type ''' not supported!']);
    end
end

if H.dims == 0
    % Empty set
    return;
end

idx = X.subs{1};
H.L = H.L(:, idx);
H.H = H.H(:, idx);

if isempty(H.L)
    H = hyperrect;
end

end
