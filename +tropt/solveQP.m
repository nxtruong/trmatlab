function [xopt, objval, success, exitflag] = solveQP(solver, Q, f, A, b, varargin)
%TROPTIM.SOLVEQP Interface to QP solvers
%   [xopt,objval,success,exitflag] = solveQP(solver,Q,f,A,b,Aeq,beq,lb,ub,x0,options)
%   Solve quadratic program (QP)
%       min_x 0.5*x'*Q*x + f'*x
%       subject to:
%           A*x <= b
%           Aeq*x = beq
%           lb <= x <= ub
%   OUTPUTS:
%     xopt   - optimal solution
%     objval - optimal value of the objective function
%     success - 0: successful, 1: infeasible, 2: unbounded, >2: other error
%     exitflag - the exit flag of the solver (see solver's docs for info)
%
%   solver is a string specifying the solver to be used; can be: quadprog,
%   gurobi.
%   If it is empty ('') then automatically chooose solver.
%
%
%See also: QUADPROG, GUROBI
%
% (C) 2012 by Truong X. Nghiem (nghiem@seas.upenn.edu)

error(nargchk(5, inf, nargin));

if isempty(solver)
    if exist('gurobi', 'file')
        solver = 'gurobi';
    elseif exist('quadprog', 'file')
        solver = 'quadprog';
    else
        error('No suitable LP solver found');
    end
end

% Process lb, ub if they are scalar
nx = length(f);
if nx > 1
    if nargin >= 8 && isscalar(varargin{3})
        % lb
        varargin{3} = repmat(varargin{3}, nx, 1);
    end
    
    if nargin >= 9 && isscalar(varargin{4})
        % ub
        varargin{4} = repmat(varargin{4}, nx, 1);
    end
end


switch solver
    case 'gurobi'
        % GUROBI
        solve_gurobi;
    case 'quadprog'
        % linprog
        solve_quadprog;
    otherwise
        error('Unknown solver for LP');
end

    function solve_quadprog
        [xopt, objval, exitflag] = quadprog(Q, f, A, b, varargin{:});
        if exitflag >= 0
            success = 0;
        elseif exitflag == -2
            success = 1;
        elseif exitflag == -3
            success = 2;
        else
            success = 3;
        end
    end

    function solve_gurobi
        % GUROBI uses column size of A for the number of variables, so even
        % if A is empty, we need to make it right.
        if isempty(A)
            A = sparse(0, length(f));
        end
        
        % Default values for optional arguments
        numvarargs = min(6, length(varargin));
        
        optargs = {[], [], [], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, LB, UB, X0, options] = optargs{:};
        
        MODEL = struct(... % MODEL
            'Q', sparse(Q/2),...   % quadratic objective matrix
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
        
        exitflag = result.status;
        switch exitflag
            case 'OPTIMAL'
                success = 0;
            case {'INFEASIBLE', 'INF_OR_UNBD'}
                success = 1;
            case 'UNBOUNDED'
                success = 2;
            otherwise
                success = 3;
        end
        
        objval = result.objval;
        xopt = result.x;
    end

end

