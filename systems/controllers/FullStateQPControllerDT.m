classdef FullStateQPControllerDT < DrakeSystem
  methods
  function obj = FullStateQPControllerDT(r,controller_data,options)
    % @param r rigid body manipulator instance
    % @param controller_data FullStateQPControllerData object containing the matrices that
    % @param options structure for specifying objective weights, slack
    % bounds, etc.
    typecheck(r,'TimeSteppingRigidBodyManipulator');
    typecheck(controller_data,'FullStateQPControllerData');
    
    if nargin>2
      typecheck(options,'struct');
    else
      options = struct();
    end
            
    input_frame = r.getStateFrame();
    output_frame = r.getInputFrame();
    
    obj = obj@DrakeSystem(0,0,input_frame.dim,output_frame.dim,true,true);
    obj = setInputFrame(obj,input_frame);
    obj = setOutputFrame(obj,output_frame);

    obj.robot = r;
    obj.numq = getNumPositions(r);
    obj.controller_data = controller_data;
    
    if isfield(options,'dt')
      % controller update rate
      typecheck(options.dt,'double');
      sizecheck(options.dt,[1 1]);
      dt = options.dt;
    else
      dt = 0.001;
    end
    obj = setSampleTime(obj,[dt;0]); % sets controller update rate
   
    % weight for grf coefficients
    if isfield(options,'w_grf')
      typecheck(options.w_grf,'double');
      sizecheck(options.w_grf,1);
      obj.w_grf = options.w_grf;
    else
      obj.w_grf = 0.0;
    end    

    % weight for cpos slack vars
    if isfield(options,'w_cpos_slack')
      typecheck(options.w_cpos_slack,'double');
      sizecheck(options.w_cpos_slack,1);
      obj.w_cpos_slack = options.w_cpos_slack;
    else
      obj.w_cpos_slack = 0.001;
    end       
    
    % weight for phi slack vars
    if isfield(options,'w_phi_slack')
      typecheck(options.w_phi_slack,'double');
      sizecheck(options.w_phi_slack,1);
      obj.w_phi_slack = options.w_phi_slack;
    else
      obj.w_phi_slack = 0.001;
    end        
    
    if isfield(options,'Kp_phi')
      typecheck(options.Kp_phi,'double');
      sizecheck(options.Kp_phi,1);
      obj.Kp_phi = options.Kp_phi;
    else
      obj.Kp_phi = 10;
    end

    if isfield(options,'Kd_phi')
      typecheck(options.Kd_phi,'double');
      sizecheck(options.Kd_phi,1);
      obj.Kd_phi = options.Kd_phi;
    else
      obj.Kd_phi = 2*sqrt(obj.Kp_phi);
    end

    % time-step for control lookahead
    if isfield(options,'timestep')
      typecheck(options.timestep,'double');
      sizecheck(options.timestep,1);
      obj.timestep = options.timestep;
    else
      obj.timestep = 0.01;
    end    
    
    
    
    % hard bound on cpos_ddot slack variables
    if isfield(options,'cpos_slack_limit')
      typecheck(options.cpos_slack_limit,'double');
      sizecheck(options.cpos_slack_limit,1);
      obj.cpos_slack_limit = options.cpos_slack_limit;
    else
      obj.cpos_slack_limit = 10;
    end

    % hard bound on phi_ddot slack variables
    if isfield(options,'phi_slack_limit')
      typecheck(options.phi_slack_limit,'double');
      sizecheck(options.phi_slack_limit,1);
      obj.phi_slack_limit = options.phi_slack_limit;
    else
      obj.phi_slack_limit = 10;
    end

    if isfield(options,'solver') 
      % 0: fastqp, fallback to gurobi barrier (default)
      % 1: gurobi primal simplex with active sets
      typecheck(options.solver,'double');
      sizecheck(options.solver,1);
      assert(options.solver==0 || options.solver==1);
    else
      options.solver = 0;
    end
    obj.solver = options.solver;
    
    if isfield(options,'offset_x')
      typecheck(options.offset_x,'logical');
      obj.offset_x = options.offset_x;
    else
      obj.offset_x = true;
    end
    
    if isfield(options,'contact_threshold')
      % minimum height above terrain for points to be in contact
      typecheck(options.contact_threshold,'double');
      sizecheck(options.contact_threshold,[1 1]);
      obj.contact_threshold = options.contact_threshold;
    else
      obj.contact_threshold = 0.001;
    end
    
    obj.gurobi_options.outputflag = 0; % not verbose
    if options.solver==0
      obj.gurobi_options.method = 2; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
    else
      obj.gurobi_options.method = 0; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
    end
    obj.gurobi_options.presolve = 0;
    % obj.gurobi_options.prepasses = 1;

    if obj.gurobi_options.method == 2
      obj.gurobi_options.bariterlimit = 20; % iteration limit
      obj.gurobi_options.barhomogeneous = 0; % 0 off, 1 on
      obj.gurobi_options.barconvtol = 5e-4;
    end
            
    [obj.jlmin, obj.jlmax] = getJointLimits(r);
  end
    
  function y=output(obj,t,~,x)
    ctrl_data = obj.controller_data;
      
    r = obj.robot;
    nq = obj.numq; 
    q = x(1:nq); 
    qd = x(nq+(1:nq)); 
            
    supp_idx = find(ctrl_data.support_times<=t,1,'last');

    h=obj.timestep;

    if ctrl_data.B_is_time_varying
      if isa(ctrl_data.B,'Trajectory')
        B_ls = fasteval(ctrl_data.B,t);
      else
        B_ls = fasteval(ctrl_data.B{supp_idx},t);
      end
    else
      B_ls = ctrl_data.B; 
    end
    R = ctrl_data.R;
    if (ctrl_data.lqr_is_time_varying)
      if isa(ctrl_data.S,'Trajectory')
        S = fasteval(ctrl_data.S,t+h);
      else
        S = fasteval(ctrl_data.S{supp_idx},t+h);
      end
      x0 = fasteval(ctrl_data.x0,t+h);
      u0 = fasteval(ctrl_data.u0,t);
    else
      S = ctrl_data.S;
      x0 = ctrl_data.x0;
      u0 = ctrl_data.u0;
    end
    q0 = x0(1:nq);

    % get phi for planned contact groups
    % if phi < threshold, then add to active contacts
    % else, add to desired contacts

    kinsol = doKinematics(r,q,true,true,qd);
    rigid_body_support_state = ctrl_data.supports(supp_idx);
    
    planned_supports = rigid_body_support_state.bodies;
    planned_contact_groups = rigid_body_support_state.contact_groups;
    %planned_num_contacts = rigid_body_support_state.num_contact_pts;      

    dim = 2; % 2D or 3D

    Jn = [];
    Jndot = [];
    phi_err = [];
    Dbar = [];
    xp = [];
    Jp = [];
    Jpdot = [];
    nc = 0;
    for j=1:length(planned_supports)
      % ridiculously inefficient for testing
      [phi,~,~,~,~,~,~,~,n,~,dn,~] = contactConstraints(r,kinsol,false,struct('terrain_only',1,...
          'body_idx',[1,planned_supports(j)],'collision_groups',planned_contact_groups(j)));
      [~,~,JB] = contactConstraintsBV(r,kinsol,false,struct('terrain_only',1,...
          'body_idx',[1,planned_supports(j)],'collision_groups',planned_contact_groups(j)));
      
      active_ind = phi<=obj.contact_threshold;
      phi_err = [phi_err;-phi(~active_ind)];
      
      nc = nc+sum(active_ind);
      Dbar = [Dbar, vertcat(JB{active_ind})']; 
      ndot = matGradMult(dn,qd);
      Jn = [Jn; n(~active_ind,:)];
      Jndot = [Jndot; ndot(~active_ind,:)];

       % hacky here because we're lacking planar system support
      terrain_pts = getTerrainContactPoints(r,planned_supports(j),planned_contact_groups(j));
      pts = [terrain_pts.pts];
      pts = pts(:,active_ind);

      xz_pts = pts([1 3],:);
      [xp_j,Jp_j] = forwardKin(r,kinsol,planned_supports(j),xz_pts,0);
      Jpdot_j = forwardJacDot(r,kinsol,planned_supports(j),pts,0);
          
      xp = [xp,xp_j];
      Jp = [Jp;Jp_j];
      Jpdot = [Jpdot;Jpdot_j];
    end
    if exist('xz_pts') && ~isempty(xz_pts)
      % compute foot placement error
      kinsol0 = r.doKinematics(q0);
      xp0 = forwardKin(r,kinsol0,planned_supports(j),xz_pts,0);
      obj.controller_data.xoffset = -1*(mean(xp0(1,:)-xp_j(1,:))); % not quite right, need to take this over all bodies in contact
    end
    
    if dim==2
       % delete y rows
      yind = 2:3:nc*3;
      Jpdot(yind,:) = [];
      Jp = sparse(Jp);
      Jpdot = sparse(Jpdot);

      nd = 2; % for friction cone approx, hard coded for now
    elseif dim==3
      nd = 4; % for friction cone approx, hard coded for now
    end
    [H,C,B] = manipulatorDynamics(r,q,qd);
    
    neps = nc*dim;
    if obj.offset_x
      xoffset = obj.controller_data.xoffset
      x0(1) = x0(1) - obj.controller_data.xoffset;
    end
    %----------------------------------------------------------------------
    % Build handy index matrices ------------------------------------------

    nu = getNumInputs(r);
    nf = nc*nd; % number of contact force variables
    neta = length(phi_err);
    nparams = nu+2*nq+nf+neps+neta;
    Iu = zeros(nu,nparams); Iu(:,1:nu) = eye(nu);
    Iq = zeros(nq,nparams); Iq(:,nu+(1:nq)) = eye(nq);
    Iqd = zeros(nq,nparams); Iqd(:,nu+nq+(1:nq)) = eye(nq);
    Ix = zeros(2*nq,nparams); Ix(:,nu+(1:2*nq)) = eye(2*nq);
    Ibeta = zeros(nf,nparams); Ibeta(:,nu+2*nq+(1:nf)) = eye(nf);
    Ieps = zeros(neps,nparams); % cpos_ddot slack vars
    Ieps(:,nu+2*nq+nf+(1:neps)) = eye(neps);
    Ieta = zeros(neta,nparams); % phi_ddot slack vars
    Ieta(:,nu+2*nq+nf+neps+(1:neta)) = eye(neta);

    
    %----------------------------------------------------------------------
    % Set up problem constraints ------------------------------------------

    lb = [r.umin' obj.jlmin' -inf*ones(1,nq) zeros(1,nf) -obj.cpos_slack_limit*ones(1,neps) -obj.phi_slack_limit*ones(1,neta)]'; % qddot/contact forces/slack vars
    ub = [r.umax' obj.jlmax' inf*ones(1,nq) inf*ones(1,nf) obj.cpos_slack_limit*ones(1,neps) obj.phi_slack_limit*ones(1,neta)]';

    Aeq_ = cell(1,4);
    beq_ = cell(1,3);

    % dynamics constraints

    h_Hinv = h*inv(H);
    if nc>0
      Aeq_{1} = Iqd - h_Hinv*B*Iu - h_Hinv*Dbar*Ibeta;
    else
      Aeq_{1} = Iqd - h_Hinv*B*Iu;
    end
    beq_{1} = qd - h_Hinv*C;

    Aeq_{2} = Iq - h*Iqd;
    beq_{2} = q;%+h*qd;

%     if nc > 0
%       % relative acceleration constraint
%       Aeq_{3} = Jp*Iqdd + Ieps;
%       beq_{3} = -Jpdot*qd - obj.Kp_accel*Jp*qd; 
%     end
% 
%     if ~isempty(phi_err)
%       phi_ddot_desired = obj.Kp_phi*phi_err - obj.Kd_phi*(Jndot*qd);
%       Aeq_{3} = Jn*Iqdd + Ieta;
%       beq_{3} = -Jndot*qd + phi_ddot_desired; 
%     end
%     
    % linear equality constraints: Aeq*alpha = beq
    Aeq = sparse(vertcat(Aeq_{:}));
    beq = vertcat(beq_{:});
    
%     phi = contactConstraints(r,x0(1:nq),false,struct('terrain_only',1));
%     active_constraints = ctrl_data.mode_data{supp_idx}.constraint_ind;
%     dynamicsFun = @(t,x,u) constrainedDynamics(r,t,x0,u,active_constraints);
%     [~,dxd] = geval(dynamicsFun,t,x0,u0,struct('grad_method','numerical'));
%     B_ls = dxd(:,getNumStates(r)+1+(1:nu));
    
    %----------------------------------------------------------------------
    % QP cost function ----------------------------------------------------
    
    Hqp = Iu'*h*R*Iu + Ix'*S*Ix;
    fqp = -x0'*S*Ix - u0'*h*R*Iu;

%     if nc > 0
%       Hqp(nu+2*nq+(1:nf),nu+2*nq+(1:nf)) = obj.w_grf*eye(nf); 
%       Hqp(nu+2*nq+nf+(1:neps),nu+2*nq+nf+(1:neps)) = obj.w_cpos_slack*eye(neps); 
%       Hqp(nu+2*nq+nf+neps+(1:neta),nu+2*nq+nf+neps+(1:neta)) = obj.w_phi_slack*eye(neta); 
%     end
%     
    %----------------------------------------------------------------------
    % Solve QP ------------------------------------------------------------

    REG = 1e-8;

    IR = eye(nparams);  
    lbind = lb>-999;  ubind = ub<999;  % 1e3 was used like inf above... right?
    Ain_fqp = full([-IR(lbind,:); IR(ubind,:)]);
    bin_fqp = [-lb(lbind); ub(ubind)];

    % call fastQPmex first
    QblkDiag = {Hqp(1:(nu+2*nq),1:(nu+2*nq)) + REG*eye(nu+2*nq), ...
                obj.w_grf*ones(nf,1) + REG*ones(nf,1), ...
                obj.w_cpos_slack*ones(neps,1) + REG*ones(neps,1), ...
                obj.w_phi_slack*ones(neta,1) + REG*ones(neta,1)};
              
    Aeq_fqp = full(Aeq);
    % NOTE: model.obj is 2* f for fastQP!!!
    [alpha,info_fqp] = fastQPmex(QblkDiag,fqp,Ain_fqp,bin_fqp,Aeq_fqp,beq,ctrl_data.qp_active_set);

    if info_fqp<0
      % then call gurobi
      disp('QPController: failed over to gurobi');
      model.Q = sparse(Hqp + REG*eye(nparams));
      model.A = Aeq;
      model.rhs = beq;
      model.sense = obj.eq_array(1:length(beq));
      model.lb = lb;
      model.ub = ub;

      model.obj = fqp;
      if obj.gurobi_options.method==2
        % see drake/algorithms/QuadraticProgram.m solveWGUROBI
        model.Q = .5*model.Q;
      end

      if (any(any(isnan(model.Q))) || any(isnan(model.obj)) || any(any(isnan(model.A))) || any(isnan(model.rhs)) || any(isnan(model.lb)) || any(isnan(model.ub)))
        keyboard;
      end

      result = gurobi(model,obj.gurobi_options);
      alpha = result.x;

      qp_active_set = find(abs(Ain_fqp*alpha - bin_fqp)<1e-6);
      obj.controller_data.qp_active_set = qp_active_set;
    end
%     beta=Ibeta*alpha
    y = Iu*alpha;
    
  end
  end

  properties (SetAccess=private)
    robot; % to be controlled
    numq;
    controller_data; % shared data handle that holds S, h, foot trajectories, etc.
    w_grf; % scalar ground reaction force weight
    w_cpos_slack; % scalar slack var weight
    w_phi_slack; % scalar slack var weight
    cpos_slack_limit; 
    phi_slack_limit; 
    Kp_phi; 
    Kd_phi; 
    gurobi_options = struct();
    solver=0;
    lc;
    eq_array = repmat('=',100,1); % so we can avoid using repmat in the loop
    ineq_array = repmat('<',100,1); % so we can avoid using repmat in the loop
    jlmin;
    jlmax;
    contact_threshold;
    offset_x; % whether or not to offset the nominal state in the x-dimension
    timestep
  end
end
