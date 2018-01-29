function v = ifelse(cond, vif, velse)
%IFELSE Functional IF-ELSE structure
%   This function returns a value according to a condition.
%     v = ifelse(cond, vif, velse)
%   If cond evaluates to true, v is assigned to vif; otherwise, it is
%   assigned to velse.
%
% (C) 2009 by Truong Nghiem (nghiem@seas.upenn.edu)

if cond
    v = vif;
else
    v = velse;
end
end