function [xopt, objval, exitflag] = solveBINLP(solver, f, A, b, varargin)
%TROPTIM.SOLVEBINLP Interface to binary LP solvers
%   [xopt,objval,exitflag] = solveBINLP(solver,f,A,b,Aeq,beq,x0,options)
%   
%   All variables are binary.  This is a special case of the MILP solved by
%   solveMILP.
%
%   OUTPUTS:
%     xopt   - optimal solution
%     objval - optimal value of the objective function
%     exitflag - 0: successful, 1: infeasible, 2: unbounded, >2: other error
%
%   solver is a string specifying the solver to be used; can be: glpk,
%   gurobi, bintprog.
%   If it is empty ('') then automatically chooose solver.
%
% (C) 2012 by Truong X. Nghiem (nghiem@seas.upenn.edu)

error(nargchk(4, inf, nargin));

if isempty(solver)
    if exist('glpk', 'file')
        solver = 'glpk';
    elseif exist('gurobi', 'file')
        solver = 'gurobi';
    elseif exist('bintprog', 'file')
        solver = 'bintprog';
    else
        error('No suitable LP solver found');
    end
end

nx = length(f);

switch solver
    case 'gurobi'
        % GUROBI
        solve_gurobi;
    case 'glpk'
        % GLPK
        solve_glpk;
    case 'bintprog'
        % bintprog
        solve_bintprog;
    otherwise
        error('Unknown solver for binary LP');
end

    function solve_glpk
        numvarargs = min(4, length(varargin));
        
        optargs = {[], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, X0, options] = optargs{:}; %#ok<ASGLU>
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
        
        [xopt, objval, exitflag] = glpk(f, A, b, [], [], ctype, repmat('B', nx, 1), 1, options);
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
        numvarargs = min(4, length(varargin));
        
        optargs = {[], [], [], []};
        optargs(1:numvarargs) = varargin(1:numvarargs);
        [Aeq, beq, X0, options] = optargs{:};
        
        MODEL = struct(... % MODEL
            'obj', f,...    % objective vector
            'A', sparse([A; Aeq]),...   % A*x <=> b
            'rhs', [b; beq],...
            'vtype', 'B');
            
        if isempty(beq)
            MODEL.sense = '<';   % no equality then use <= for all
        else
            MODEL.sense = [repmat('<', 1, length(b)), repmat('=', 1, length(beq))];
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

    function solve_bintprog
        [xopt, objval, exitflag] = bintprog(f, A, b, varargin{:});
        if exitflag >= 0
            exitflag = 0;
        elseif exitflag == -2
            exitflag = 1;
        else
            exitflag = 3;
        end
    end
end

