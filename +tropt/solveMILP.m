function [xopt, objval, exitflag] = solveMILP(solver, f, A, b, varargin)
%TROPTIM.SOLVEMILP Interface to MILP solvers
%   [xopt,objval,exitflag] = solveMILP(solver,f,A,b,Aeq,beq,lb,ub,vartype,x0,options)
%   
%   vartype is a string of 'C' or 'B' or 'I'
%
%   OUTPUTS:
%     xopt   - optimal solution
%     objval - optimal value of the objective function
%     exitflag - 0: successful, 1: infeasible, 2: unbounded, >2: other error
%
%   solver is a string specifying the solver to be used; can be: glpk,
%   gurobi.
%   If it is empty ('') then automatically chooose solver.
%
% (C) 2012 by Truong X. Nghiem (nghiem@seas.upenn.edu)

error(nargchk(4, inf, nargin));

if isempty(solver)
    if exist('gurobi', 'file')
        solver = 'gurobi';
    elseif exist('glpk', 'file')
        solver = 'glpk';
    else
        error('No suitable LP solver found');
    end
end

% Process lb, ub if they are scalar
nx = length(f);
if nx > 1
    if nargin >= 7 && isscalar(varargin{3})
        % lb
        varargin{3} = repmat(varargin{3}, nx, 1);
    end
    
    if nargin >= 8 && isscalar(varargin{4})
        % ub
        varargin{4} = repmat(varargin{4}, nx, 1);
    end
end


switch solver
    case 'gurobi'
        % GUROBI
        solve_gurobi;
    case 'glpk'
        % GLPK
        solve_glpk;
    otherwise
        error('Unknown solver for MILP');
end

    function solve_glpk
        numvarargs = min(7, length(varargin));
        
        optargs = {[], [], [], [], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, LB, UB, vartype, X0, options] = optargs{:}; %#ok<ASGLU>
        % X0 is not used by GLPK
        
        if ~isempty(Aeq) && ~isempty(beq)
            A = [Aeq; A];
            b = [beq; b];
            ctype = blanks(length(b));
            leq = length(beq);
            ctype(1:leq) = 'S';
            ctype(leq+1:end) = 'U';
        else
            ctype = [];
        end
        
        if isempty(options)
            options = struct('msglev', 0);  % turn off printing out
        elseif ~isfield(options, 'msglev')
            options.msglev = 0;
        end
        
        if isscalar(vartype)
            vartype = repmat(vartype, nx, 1);
        else
            vartype = vartype(:);
        end
        
        [xopt, objval, exitflag] = glpk(f, A, b, LB, UB, ctype, vartype, 1, options);
        switch exitflag
            case 5
                exitflag = 0;
            case 4
                exitflag = 1;
            case 6
                exitflag = 2;
            otherwise
                exitflag = 3;
        end
    end

    function solve_gurobi
        % GUROBI uses column size of A for the number of variables, so even
        % if A is empty, we need to make it right.
        if isempty(A)
            A = sparse(0, length(f));
        end
        
        % Default values for optional arguments
        numvarargs = min(7, length(varargin));
        
        optargs = {[], [], [], [], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, LB, UB, vartype, X0, options] = optargs{:};
        
        MODEL = struct(... % MODEL
            'obj', f,...    % objective vector
            'A', sparse([A; Aeq]),...   % A*x <=> b
            'rhs', [b; beq]);
            
        if isempty(beq)
            MODEL.sense = '<';   % no equality then use <= for all
        else
            MODEL.sense = [repmat('<', 1, length(b)), repmat('=', 1, length(beq))];
        end
        
        if isempty(LB)
            % GUROBI by default sets LB = 0 if no bounds are given
            % What we want is -inf
            MODEL.lb = -inf(nx, 1);
        else
            MODEL.lb = LB;
        end
        
        % By default, GUROBI uses upper bound of +inf, so we don't need to
        % set
        if ~isempty(UB)
            MODEL.ub = UB;
        end
        
        % vartype
        if ~isempty(vartype)
            MODEL.vtype = vartype(:);
        end
        
        % X0
        if ~isempty(X0)
            MODEL.start = double(X0);
        end
        
        if isempty(options)
            % Shut off GUROBI
            options = struct('OutputFlag', 0);
        elseif ~isfield(options, 'OutputFlag')
            options.OutputFlag = 0;
        end
        
        result = gurobi(MODEL, options);
        
        switch result.status
            case 'OPTIMAL'
                exitflag = 0;
            case {'INFEASIBLE', 'INF_OR_UNBD'}
                exitflag = 1;
            case 'UNBOUNDED'
                exitflag = 2;
            otherwise
                exitflag = 3;
        end
        
        objval = result.objval;
        xopt = result.x;
    end
end

