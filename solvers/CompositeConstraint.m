classdef CompositeConstraint
  %CompositeConstraint
  % This is a container class for a set of general constraints
  % and slack variables.
  % The idea is that some intuitive constraints are best formatted as a
  % related set of constraints (and, that there may be multiple ways to
  % describe these, which the user may wish to switch back and forth
  % between)
  
  properties
    constraints = {}
    n_slack = 0
    slack_name
    costs = {}
  end
  
  methods
    function obj = CompositeConstraint(constraints,n,slack_name,costs)
      if nargin > 0
        if ~iscell(constraints)
          if ~isa(constraints,'Constraint')
            error('Drake:CompositeConstraint:InvalidArgument','constraints argument must be a cell-array of constraints or a Constraint type');
          end
          constraints = {constraints};
        end
        obj.constraints = constraints;
      end
      if nargin > 1
        obj.n_slack = n;
      end
      if nargin > 2
        if(~iscellstr(slack_name))
          error('Drake:CompositeConstraint:InvalidArgument','slack_name argument must be a cell of strings');
        end
        if(numel(slack_name) ~= n)
          error('Drake:CompositeConstraint:InvalidArgument','slack_name should have %d elements',obj.n_slack);
        end
        slack_name = slack_name(:);
      else
        slack_name = repmat({'slack'},obj.n_slack,1);
      end
      obj.slack_name = slack_name;
      if nargin > 3
        if ~iscell(costs)
          if ~isa(costs,'Constraint')
            error('Drake:CompositeConstraint:InvalidArgument','costs argument must be a cell array of constraints or a Constraint type');
          end
          costs = {costs};
        end
        obj.costs = costs;
      end
    end

    
    function obj = addConstraints(obj, constraints)
      obj.constraints = [obj.constraints;constraints];
    end
    
    function obj = addSlackVariables(obj, n,slack_name)
      obj.n_slack = obj.n_slack + n;
      if(nargin>2)
        if(~iscellstr(slack_name) || numel(slack_name) ~= n)
          error('Drake:CompositeConstraint:InvalidArgument','slack_name should be a cell of strings, with number of elements equals %d',n);
        end
        slack_name = slack_name(:);
      else
        slack_name = repmat({'slack'},n);
      end
      obj.slack_name = [obj.slack_name;slack_name];
    end
    
    function obj = addCosts(obj, costs)
      obj.costs = [obj.costs; costs];
    end
  end
  
end

 

 