function [xopt, objval, exitflag] = solveLP(solver, f, A, b, varargin)
%TROPTIM.SOLVELP Interface to LP solvers
%   [xopt,objval,exitflag] = solveLP(solver,f,A,b,Aeq,beq,lb,ub,options)
%   See linprog for input arguments.
%   OUTPUTS:
%     xopt   - optimal solution
%     objval - optimal value of the objective function
%     exitflag - 0: successful, 1: infeasible, 2: unbounded, >2: other error
%
%   solver is a string specifying the solver to be used; can be: cdd, glpk,
%   linprog, clp, gurobi.
%   If it is empty ('') then automatically chooose solver.
%
%   For cdd, if options = 'DS' then the Dual-Simplex method is used.
%
%See also: LINPROG
%
% (C) 2018 by Truong X. Nghiem (nghiem@seas.upenn.edu)

% History:
% 2018-01-29 Truong removed X0 argument; it's not used by any solver.

narginchk(4, inf);

if isempty(solver)
    if exist('cddmex', 'file')
        solver = 'cdd';
    elseif exist('gurobi', 'file')
        solver = 'gurobi';
    elseif exist('glpk', 'file')
        solver = 'glpk';
    elseif exist('clp', 'file')
        solver = 'clp';
    elseif exist('linprog', 'file')
        solver = 'linprog';
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
    case 'cdd'
        % CDD
        solve_cdd;
    case 'gurobi'
        % GUROBI
        solve_gurobi;
    case 'glpk'
        % GLPK
        solve_glpk;
    case 'clp'
        % CLP
        solve_clp;
    case 'linprog'
        % linprog
        solve_linprog;
    otherwise
        error('Unknown solver for LP');
end

    function solve_linprog
        [xopt, objval, exitflag] = linprog(f, A, b, varargin{:});
        if exitflag >= 0
            exitflag = 0;
        elseif exitflag == -2
            exitflag = 1;
        elseif exitflag == -3
            exitflag = 2;
        else
            exitflag = 3;
        end
    end

    function solve_cdd
        ndim = length(f);
        numvarargs = min(5, length(varargin));
        
        % Default values for optional arguments
        optargs = {[], [], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, LB, UB, options] = optargs{:}; %#ok<NASGU>
        
        IN = struct('obj', f.', 'A', [Aeq; A], 'B', [beq; b], 'lin', 1:length(beq));
        if ~isempty(LB)
            IN.A = [IN.A; -eye(ndim)];
            IN.B = [IN.B; -LB];
        end
        if ~isempty(UB)
            IN.A = [IN.A; eye(ndim)];
            IN.B = [IN.B; UB];
        end

        if isempty(options)
            OUT = cddmex('solve_lp', IN);
        else
            OUT = cddmex('solve_lp_DS', IN);
        end
        xopt = OUT.xopt;
        objval = OUT.objlp;
        if OUT.how == 1
            exitflag = 0;
        elseif OUT.how == 6
            exitflag = 2;
        elseif OUT.how == 2
            exitflag = 1;
        else
            exitflag = 3;
        end
    end

    function solve_glpk
        numvarargs = min(5, length(varargin));
        
        % Default values for optional arguments
        optargs = {[], [], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, LB, UB, options] = optargs{:}; %#ok<ASGLU>
        
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
        
        [xopt, objval, exitflag] = glpk(f, A, b, LB, UB, ctype, [], 1, options);
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
        numvarargs = min(5, length(varargin));
        
        optargs = {[], [], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, LB, UB, options] = optargs{:};
        
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
        
        if isfield(result, 'objval')
            objval = result.objval;
        else
            objval = [];
        end
        if isfield(result, 'x')
            xopt = result.x;
        else
            xopt = [];
        end
    end
    
    function solve_clp
        numvarargs = length(varargin);
        if numvarargs > 4, numvarargs = 4; end  % Only use the first 4
        [xopt, ~, exitflag] = clp([], f, A, b, varargin{1:numvarargs});
        objval = f'*xopt;
    end

end

