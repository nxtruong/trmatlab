function [Ain,bin,Aeq,beq] = getconstraints(obj)
%GETCONSTRAINTS Get linear constraints of the given hyperrectangle.
%   Returns the inequality and equality constraints describing the given
%   hyperrectangle.
%
%       [Ain,bin,Aeq,beq] = getconstraints(H)
%
%   H must be a single hyperrectangle.
%   H = {x | Ain*x <= bin and Aeq*x = beq}
%
% (C) 2011 by Truong X. Nghiem (nghiem@seas.upenn.edu)

assert(length(obj) == 1, 'A single hyperrectangle must be provided.');

% Initialize outputs
if obj.dims == 0
    Ain = [];
    bin = [];
    Aeq = [];
    beq = [];
    return;
end

% Find the number of inequalities and equalities
idxeq = obj.L == obj.H;
neq = sum(idxeq);
idxin = ~idxeq;

% Construct inequalities
if neq < obj.dims
    Ain = eye(obj.dims);
    Ain(idxeq,:) = [];    % Remove dimensions c.t. equalities
    Ain = [Ain; -Ain];
    bin = [obj.H(idxin); -obj.L(idxin)];
else
    Ain = [];
    bin = [];
end

if neq > 0
    Aeq = eye(obj.dims);
    Aeq(idxin,:) = [];
    beq = obj.L(idxeq);
else
    Aeq = [];
    beq = [];
end

end
