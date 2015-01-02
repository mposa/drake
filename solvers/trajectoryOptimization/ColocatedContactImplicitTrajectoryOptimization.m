classdef ColocatedContactImplicitTrajectoryOptimization < DirectTrajectoryOptimization
  % phi, lambda
  properties
    nC
    nD % number of friction elements per contact
    
    lp_inds % post-t orderered [lambda_N;lambda_f1;lambda_f2;...] for each contact sequentially
    lm_inds % orderered [lambda_N;lambda_f1;lambda_f2;...;] for each contact sequentially
    L_inds % orderered [Lambda_N;Lambda_f1;Lambda_f2;...;] for each contact sequentially
    gammap_inds;
    gammam_inds;
    dgamma_inds;
    %     Gamma_inds;
    phi_inds % gap function, ordered sequentially
    eta_inds % post-impact normal velocity
    
    ljlp_inds  % joint limit forces
    ljlm_inds  % joint limit forces
    Ljl_inds  % joint limit forces
    etajl_inds  % joint limit forces
    
    jl_lb_ind  % joint indices where the lower bound is finite
    jl_ub_ind % joint indices where the lower bound is finite
    nJL % number of joint limits = length([jl_lb_ind;jl_ub_ind])
  end
  
  
  methods
    function obj = ColocatedContactImplicitTrajectoryOptimization(plant,N,duration,options)
      if nargin<4, options=struct(); end
      
      if ~isfield(options,'nlcc_mode')
        options.nlcc_mode = 2;
      end
      if ~isfield(options,'lincc_mode')
        options.lincc_mode = 1;
      end
      if ~isfield(options,'compl_slack')
        options.compl_slack = 0;
      end
      if ~isfield(options,'lincompl_slack')
        options.lincompl_slack = 0;
      end
      if ~isfield(options,'jlcompl_slack')
        options.jlcompl_slack = 0;
      end
      if ~isfield(options,'lambda_mult')
        options.lambda_mult = 1;
      end
      if ~isfield(options,'Lambda_mult')
        options.Lambda_mult = 1;
      end
      if ~isfield(options,'lambda_jl_mult')
        options.lambda_jl_mult = 1;
      end
      if ~isfield(options,'active_collision_options')
        options.active_collision_options.terrain_only = true;
      end
      
      obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
      
    end
    
    function obj = addDynamicConstraints(obj)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      nq = obj.plant.getNumPositions();
      N = obj.N;
      
      constraints = cell(N-1,1);
      lincompl_constraints = cell(N-1,1);
      nonlincompl_constraints = cell(N-1,1);
      jlcompl_constraints = cell(N-1,1);
      philcp_constraints = cell(N,1);
      Lzlcp_constraints = cell(N,1);
      dyn_inds = cell(N-1,1);
      
      nDynVars = 1 + 2*nX + 2*nU + 3*obj.nC*(1+obj.nD) + 2*obj.nJL;
      
      cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),nDynVars,@obj.dynamics_wimpulse_constraint_fun);
      phi_cnstr = FunctionHandleConstraint(zeros(obj.nC,1),zeros(obj.nC,1),nq + obj.nC,@obj.phi_fun);
      eta_cnstr =  FunctionHandleConstraint(zeros(2*obj.nC,1),inf(2*obj.nC,1),nX + obj.nC*(2+obj.nD),@obj.eta_fun);
      gamma_cnstr =  FunctionHandleConstraint(zeros(2*obj.nC*obj.nD,1),inf(2*obj.nC*obj.nD,1),nX+obj.nC*(3+obj.nD),@obj.gamma_fun);
      dgamma_cnstr =  FunctionHandleConstraint(zeros(obj.nD*obj.nC,1),inf(obj.nD*obj.nC,1),2*nq+obj.nC,@obj.dgamma_fun);
      
      nGammacVars = 1 + 2*nX + 2*nU + 3*obj.nC*(1+obj.nD) + 2*obj.nJL + obj.nC;
      gammac_cnstr = FunctionHandleConstraint(zeros(obj.nD*obj.nC,1),inf(obj.nD*obj.nC,1),nGammacVars,@obj.gammac_fun);
      
      [~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(zeros(nq,1),false,obj.options.active_collision_options);
      mu = 1;
      
      %       args = {rand;randn(12,1);randn(12,1);randn(3,1);randn(3,1);randn(6,1);zeros(0,1);randn(6,1);zeros(0,1);randn(6,1)};
      % x0 = zeros(12,1);
      %       x0(2) = 1;
      % args = {3;x0+1e-4*randn(12,1);x0;zeros(3,1);zeros(3,1);[.5;0;0;.5;0;0];zeros(0,1);[.5;0;0;.5;0;0];zeros(0,1);zeros(6,1)};
      %       [g,dg] = dynamics_constraint_fun(obj,args{1:end-1});
      %       [g,dg2] = geval(@obj.dynamics_constraint_fun,args{1:end-1},struct('grad_method','numerical'));
      %
      
      % friction limits
      A_friction = kron(eye(obj.nC),[1 -ones(1,obj.nD)]);
      A_friction(:,1:1+obj.nD:end) = diag(mu*ones(obj.nC,1));
      fric_cnstr = LinearConstraint(zeros(obj.nC,1),inf(obj.nC,1),A_friction);
      
      % shared data functions
      shared_data_index = obj.getNumSharedDataFunctions;
      for i=1:obj.N,
        obj = obj.addSharedDataFunction(@obj.contact_data_fun,{obj.x_inds(1:nq,i);obj.x_inds(nq+1:end,i);obj.L_inds(:,i)});
      end
      
      for i=1:obj.N-1,
        %         dyn_inds{i} = [obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.l_inds(:,i);obj.ljl_inds(:,i)];
        dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1);...
          obj.lp_inds(:,i);obj.ljlp_inds(:,i);obj.lm_inds(:,i+1);obj.ljlm_inds(:,i+1);obj.L_inds(:,i)};
        constraints{i} = cnstr;
        
        for j=1:length(dyn_inds{i}),
          args{j} = randn(length(dyn_inds{i}{j}),1);
        end
        
        obj = obj.addConstraint(constraints{i}, dyn_inds{i}, [shared_data_index+i; shared_data_index+i+1]);
        
      end
      
      
      for i=1:obj.N,
        if obj.nC > 0
          % bounding box constraints on forces
          bbox_inds = [obj.lp_inds(:);obj.lm_inds(:);obj.L_inds(:);obj.phi_inds(:);obj.gammap_inds(:);obj.gammam_inds(:);obj.dgamma_inds(:)];
          nbbox = length(bbox_inds);
          positivity_cnstr = BoundingBoxConstraint(zeros(nbbox,1),inf(nbbox,1));
          obj = obj.addConstraint(positivity_cnstr,bbox_inds);
          
          %phi slack variable
          obj = obj.addConstraint(phi_cnstr,{obj.x_inds(1:nq,i);obj.phi_inds(:,i)},shared_data_index+i);
          
          % eta slack variable
          obj = obj.addConstraint(eta_cnstr,{obj.x_inds(1:nq,i);obj.x_inds(nq+1:end,i);obj.L_inds(:,i);obj.eta_inds(:,i)},shared_data_index+i);
          
          % gamma slack variables
          obj = obj.addConstraint(gamma_cnstr,{obj.x_inds(1:nq,i);obj.x_inds(nq+1:end,i);obj.L_inds(:,i);obj.gammam_inds(:,i);obj.gammap_inds(:,i)},shared_data_index+i);
          %               function [f,df] = dgamma_fun(obj,q0,q1,dgamma,data0,data1)
          
          
          
          use_gammac = 1;
          %using dgamma OR gammac
          if use_gammac
            if i~=obj.N
              gammac_cnstr_inds = {dyn_inds{i}{:},obj.dgamma_inds(:,i)};
              obj = obj.addConstraint(gammac_cnstr,gammac_cnstr_inds,[shared_data_index+i;shared_data_index+i+1]);
            end
          else
            if i~=obj.N
              %dgamma slack
              %gammac_fun(obj,h,x0,x1,u0,u1,lambda0,lambda_jl0,lambda1,lambda_jl1,gammac,data0,data1)
              obj = obj.addConstraint(dgamma_cnstr,{obj.x_inds(1:nq,i);obj.x_inds(1:nq,i+1);obj.dgamma_inds(:,i)},[shared_data_index+i;shared_data_index+i+1]);
            end
          end
          
          
          % friction limits
          obj = obj.addConstraint(fric_cnstr,obj.lm_inds(:,i));
          obj = obj.addConstraint(fric_cnstr,obj.lp_inds(:,i));
          obj = obj.addConstraint(fric_cnstr,obj.L_inds(:,i));
          
          
          % Linear complementarity constraints
          % phi_k /perp lambdap_k,z + lambdam_k,z + lambdap_k-1,z + lambdam_k+1,z
          lz_inds = 1:1+obj.nD:size(obj.L_inds,1);
          
          philcp_vars = [obj.lp_inds(lz_inds,i);obj.lm_inds(lz_inds,i);obj.phi_inds(:,i)];
          if i==1,
            M_phi = repmat(eye(obj.nC),1,3);
            philcp_vars = [obj.lm_inds(lz_inds,i+1);philcp_vars];
          elseif i==obj.N
            M_phi = repmat(eye(obj.nC),1,3);
            philcp_vars = [obj.lp_inds(lz_inds,i-1);philcp_vars];
          else
            M_phi = repmat(eye(obj.nC),1,4);
            philcp_vars = [obj.lp_inds(lz_inds,i-1);obj.lm_inds(lz_inds,i+1);philcp_vars];
          end
          r_phi = zeros(obj.nC,1);
          W_phi = zeros(obj.nC);
          philcp_constraints{i} = LinearComplementarityConstraint(W_phi,r_phi,M_phi,obj.options.lincc_mode,obj.options.lincompl_slack);
          obj = obj.addConstraint(philcp_constraints{i},philcp_vars);
          
          % Lambda_k,z /perp phi_k + eta_k
          r_Lz = zeros(obj.nC,1);
          W_Lz = zeros(obj.nC);
          M_Lz = [eye(obj.nC) eye(obj.nC)];
          lzlcp_vars = [obj.phi_inds(:,i);obj.eta_inds(:,i);obj.L_inds(lz_inds,i)];
          Lzlcp_constraints{i} = LinearComplementarityConstraint(W_Lz,r_Lz,M_Lz,obj.options.lincc_mode,obj.options.lincompl_slack);
          obj = obj.addConstraint(Lzlcp_constraints{i},lzlcp_vars);
          
          
          %noslip constraint
          %gamma /perp lambda_z
%           W_ns = zeros(obj.nC);
%           r_ns = zeros(obj.nC,1);
%           M_ns = [eye(obj.nC)];
%           ns_constraint = LinearComplementarityConstraint(W_ns,r_ns,M_ns,obj.options.lincc_mode,obj.options.lincompl_slack);
%           obj = obj.addConstraint(ns_constraint,[obj.L_inds(lz_inds,i);obj.gammap_inds(:,i)]);
          M_ns = [eye(obj.nC) zeros(obj.nC) zeros(obj.nC); zeros(obj.nC) eye(obj.nC) eye(obj.nC);];
          r_ns = zeros(2*obj.nC,1);
          W_ns = zeros(2*obj.nC);
          ns_vars = [obj.lm_inds(lz_inds,i);obj.lp_inds(lz_inds,i);obj.L_inds(lz_inds,i)];
          if i~= obj.N,
            M_ns = [M_ns [zeros(obj.nC);eye(obj.nC)]];
            ns_vars = [ns_vars;obj.lm_inds(lz_inds,i+1)];
          end
          
          if i~= 1
            M_ns = [M_ns [eye(obj.nC);zeros(obj.nC)]];
            ns_vars = [ns_vars;obj.lp_inds(lz_inds,i-1)];
          end
          ns_vars = [ns_vars;obj.gammam_inds(:,i);obj.gammap_inds(:,i)];
          
          ns_constraint = LinearComplementarityConstraint(W_ns,r_ns,M_ns,obj.options.lincc_mode,obj.options.lincompl_slack);
          obj = obj.addConstraint(ns_constraint,ns_vars);
          
          if i~= obj.N
            %dgamma constraint
            %dgamma /perp lambdazp + lambdazm
            W_dg = zeros(obj.nC);
            r_dg = zeros(obj.nC,1);
            M_dg = [eye(obj.nC) eye(obj.nC)];
            dg_constraint = LinearComplementarityConstraint(W_dg,r_dg,M_dg,obj.options.lincc_mode,obj.options.lincompl_slack);
            obj = obj.addConstraint(dg_constraint,[obj.lp_inds(lz_inds,i);obj.lm_inds(lz_inds,i+1);obj.dgamma_inds(:,i)]);
          end
        end
        
        if obj.nJL > 0
          % joint limit linear complementarity constraint
          % lambda_jl /perp [q - lb_jl; -q + ub_jl]
          W_jl = zeros(obj.nJL);
          [r_jl,M_jl] = jointLimitConstraints(obj.plant,zeros(nq,1));
          jlcompl_constraints{i} = LinearComplementarityConstraint(W_jl,r_jl,M_jl,obj.options.lincc_mode,obj.options.jlcompl_slack);
          
          obj = obj.addConstraint(jlcompl_constraints{i},[obj.x_inds(1:nq,i+1);obj.ljl_inds(:,i)]);
        end
      end
    end
    
    function [f,df] = dgamma_fun(obj,q0,q1,dgamma,data0,data1)
      % dgamma - tang(0) + tang(1)
      
      f = kron(dgamma,ones(obj.nD,1)) - data0.tang_pos(:) + data1.tang_pos(:);
      df = [-data0.dtang_pos data1.dtang_pos kron(eye(length(dgamma)),ones(obj.nD,1))];      
    end
    
    function [f,df] = gamma_fun(obj,q,vm,Lambda,gammam,gammap,data)
      % f= [gammam + D*psim;gammap + D*psip]
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nl = length(Lambda);
      use_data = (nargin > 6);
      fm = zeros(obj.nC*obj.nD,1);
      dfm = zeros(obj.nC*obj.nD,nq+nv+obj.nC*(3+obj.nD));
      fp = zeros(obj.nC*obj.nD,1);
      dfp = zeros(obj.nC*obj.nD,nq+nv+obj.nC*(3+obj.nD));
      
      if use_data
        vp = data.vp;
        dvp = data.dvp;
        D = data.D;
        dD = data.dD;
      else
        error('Not implemented yet, need to copy/paste')
      end
      
      for j=1:obj.nD,
        fm(j:obj.nD:end) = gammam+D{j}*vm;
        dfm(j:obj.nD:end,nq+nv+nl+(1:obj.nC)) = eye(size(D{j},1));  %d/dgammam
        dfm(j:obj.nD:end,nq+(1:nv)) = D{j};%d/dv
        dfm(j:obj.nD:end,1:nq) = matGradMult(dD{j},vm);%d/dq
        
        fp(j:obj.nD:end) = gammap+D{j}*vp;
        dfp(j:obj.nD:end,nq+nv+nl+obj.nC+(1:obj.nC)) = eye(size(D{j},1));  %d/dgammap
        
        dfp(j:obj.nD:end,1:nq+nv+nl) = D{j}*dvp;
        dfp(j:obj.nD:end,1:nq) = dfp(j:obj.nD:end,1:nq)+ matGradMult(dD{j},vp);%d/dq
      end
      
      f = [fm;fp];
      df = [dfm;dfp];
    end
    
    function [f,df] = gammac_fun(obj,h,x0,x1,u0,u1,lambda0,lambda_jl0,lambda1,lambda_jl1,Lambda,gammac,data0,data1)
      use_data = (nargin > 11);
      
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      nl = length(lambda0);
      njl = length(lambda_jl0);
      
      q0 = x0(1:nq);
      
      if use_data
        v0p = data0.vp;
        dv0pdLambda = data0.dvp(:,nX+1:end);
        dv0pdx0 = data0.dvp(:,1:nX);
        dx0pdx0 = [eye(nq) zeros(nq,nv); dv0pdx0];
        x0p = [q0;v0p];
        
        [xdot0,dxdot0p] = constrained_dynamics(obj,x0p(1:nq),x0p(nq+1:end),u0,lambda0,lambda_jl0,data0);
        v0p_inds = 1+nq:nq+nv;
        dxdot0 = [dxdot0p dxdot0p(:,v0p_inds)*dv0pdLambda];
        dxdot0(:,2:1+nX) = dxdot0(:,2:1+nX)*dx0pdx0;
        
        %dxdot0 and dxdot1 different sizes, clearly something's wrong here
        
        [xdot1,dxdot1] = constrained_dynamics(obj,x1(1:nq),x1(nq+1:end),u1,lambda1,lambda_jl1,data1);
      else
        error('Not implemented')               
      end
                  
      
      % cubic interpolation to get xcol and xdotcol, as well as
      % derivatives
      xcol = .5*(x0p+x1) + h/8*(xdot0-xdot1);
      dxcol = [1/8*(xdot0-xdot1), (.5*dx0pdx0 + h/8*dxdot0(:,1:nX)), ...
        (.5*eye(nX) - h/8*dxdot1(:,1:nX)), h/8*dxdot0(:,nX+1:nX+nU) -h/8*dxdot1(:,nX+1:nX+nU), ...
        h/8*dxdot0(:,nX+nU+1:nX+nU+nl+njl), -h/8*dxdot1(:,nX+nU+1:end), h/8*dxdot0(:,nX+nU+nl+njl+1:end)+.5*[zeros(nq,nl); dv0pdLambda]];
      
%       dxcol(:,2+2*nX+2*nU+2*nl+2*njl:end) = dxcol(:,2+2*nX+2*nU+2*nl+2*njl:end) + dxcol(:,2+nq:1+nX)*dv0pdLambda;
%       dxcol(:,2:1+nX) = dxcol(:,2:1+nX)*dx0pdx0;
      
      qc = xcol(1:nq);
      vc = xcol(nq+1:end);
      
      [phi,~,d,xA,xB,idxA,idxB,~,n,D,dn,dD] = obj.plant.contactConstraints(qc,false,obj.options.active_collision_options);
      
      f = zeros(obj.nC*obj.nD,1);
      df = zeros(obj.nC*obj.nD,1+2*nX+2*nU+3*nl+2*njl+obj.nC);
      
      for j=1:obj.nD,
        f(j:obj.nD:end) = gammac + D{j}*vc;
        
        df(j:obj.nD:end,end-obj.nC+1:end) = eye(size(D{j},1));  %d/dgammac                
        df(j:obj.nD:end,1:end-obj.nC) = D{j}*dxcol(nq+1:end,:) + matGradMult(dD{j}*dxcol(1:nq,:),vc);%d/dv
      end      
    end
    
    function data = contact_data_fun(obj,q,vm,Lambda)
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nl = length(Lambda);
      
      [H,C,B,dH] = obj.plant.manipulatorDynamics(q,vm);
      Hinv = inv(H)';
      [phi,~,d,xA,xB,idxA,idxB,~,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
      J = zeros(nl,nq);
      dJ = zeros(nl*nq,nq);
      
      if obj.nC > 0
        J(1:1+obj.nD:end,:) = n;
        
        dJ(1:1+obj.nD:end,:) = dn;
        
        for j=1:length(D),
          J(1+j:1+obj.nD:end,:) = D{j};
          dJ(1+j:1+obj.nD:end,:) = dD{j};
        end
        
        vp = vm + Hinv*J'*Lambda;
        dvp = [-Hinv*matGradMult(dH(:,1:nq),Hinv*J'*Lambda)+Hinv*matGradMult(dJ,Lambda,true), eye(nv), Hinv*J'];
      else
        vp = vm;
        dvp = [zeros(nv,nq) eye(nv)];
      end      
      
      
      data.dvp = dvp;
      data.vp = vp;
      data.Hinv = Hinv;
      data.H = H;
      data.dH = dH;
      data.phi = phi;
      data.n = n;
      data.dn = dn;
      data.dD = dD;
      data.D = D;
      data.J = J;
      data.dJ = dJ;
      
      body_inds = [idxA(idxA ~= 1);idxB(idxB ~= 1)];
      xBody = [xA(:,idxA ~= 1) xB(idxB ~= 1)];
      kinsol = obj.plant.doKinematics(q);      
      
      tang_pos = zeros(2*length(d),size(xBody,2));
      dtang_pos = zeros(2*length(d)*size(xBody,2),nq);
      for i=1:size(xBody,2),
        [xWorld, dxWorld] = obj.plant.forwardKin(kinsol,body_inds(i),xBody(:,i));
        for j=1:length(d),
          tang_pos(2*j-1,i) = xWorld'*d{j}(:,i);
          tang_pos(2*j,i) = -xWorld'*d{j}(:,i);
          dtang_pos(2*j - 1 + 2*(i-1)*length(d),:) = d{j}(:,i)'*dxWorld;
          dtang_pos(2*j + 2*(i-1)*length(d),:) = -d{j}(:,i)'*dxWorld;
        end
      end

      data.tang_pos = tang_pos;
      data.dtang_pos = dtang_pos;      
    end
    
    function [f,df] = eta_fun(obj,q,vm,Lambda,eta,data)
      nEta = length(eta);
      nv = obj.plant.getNumVelocities;
      nl = length(Lambda);
      use_data = (nargin > 5);
      
      if use_data
        n = data.n;
        dn = data.dn;
        vp = data.vp;
        dvp = data.dvp;
      else
        % just need to copy in these calculations
        error('Not implemented yet')
      end
      
      phidot = n*vp;
      dphidot = [matGradMult(dn,vp) zeros(nEta,nv+nl)] + n*dvp;
      
      f = [eta - phidot; eta + phidot];
      df = [-dphidot eye(nEta);dphidot eye(nEta)];
    end
    
    function [f,df] = phi_fun(obj,q,phi_var,data)
      use_data = (nargin > 3);
      
      if use_data
        phi = data.phi;
        dphi = data.n;
      else
        [phi,~,~,~,~,~,~,~,dphi] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
      end
      f = phi - phi_var;
      df = [dphi -eye(length(phi))];
    end
    
    function [f,df] = dynamics_wimpulse_constraint_fun(obj,h,x0,x1,u0,u1,lambda0,lambda_jl0,lambda1,lambda_jl1,Lambda,data0,data1)
      use_data = (nargin > 11);
      
      lambda0 = lambda0*obj.options.lambda_mult;
      lambda1 = lambda1*obj.options.lambda_mult;
      %       Lambda = Lambda*obj.options.Lambda_mult;
      
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      nl = length(lambda0);
      njl = length(lambda_jl0);
      q0 = x0(1:nq);
      v0m = x0(nq+1:end);
      
      % calculate v0p, velocity after impact
      if use_data
        v0p = data0.vp;
        dv0pdLambda = data0.dvp(:,nX+1:end);
        dv0pdx0 = data0.dvp(:,1:nX);
      else
        [H,C,B,dH] = obj.plant.manipulatorDynamics(q0,v0m);
        [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q0,false,obj.options.active_collision_options);
        J = zeros(nl,nq);
        J(1:1+obj.nD:end,:) = n;
        dJ = zeros(nl*nq,nq);
        dJ(1:1+obj.nD:end,:) = dn;
        for j=1:length(D),
          J(1+j:1+obj.nD:end,:) = D{j};
          dJ(1+j:1+obj.nD:end,:) = dD{j};
        end
        Hinv = inv(H);
        v0p = v0m + Hinv*J'*Lambda;
        dv0pdLambda = Hinv*J';
        dv0pdx0 = [-Hinv*matGradMult(dH(:,1:nq),Hinv*J'*Lambda)+Hinv*matGradMult(dJ,Lambda,true), eye(nv)];
      end
      
      dx0pdx0 = [eye(nq) zeros(nq,nv); dv0pdx0];
      x0p = [q0;v0p];
      
      if use_data
        [g,dg] = dynamics_constraint_fun(obj,h,x0p,x1,u0,u1,lambda0,lambda_jl0,lambda1,lambda_jl1,data0,data1);
      else
        [g,dg] = dynamics_constraint_fun(obj,h,x0p,x1,u0,u1,lambda0,lambda_jl0,lambda1,lambda_jl1);
      end
      
      v0p_inds = 2+nq:1+nq+nv;
      
      f = g;
      df = [dg dg(:,v0p_inds)*dv0pdLambda];
      df(:,2:1+nX) = df(:,2:1+nX)*dx0pdx0;
      
      df(:,2+2*nX+2*nU:1+2*nX+2*nU+nl) = df(:,2+2*nX+2*nU:1+2*nX+2*nU+nl)*obj.options.lambda_mult;
      df(:,2+2*nX+2*nU+nl+njl:1+2*nX+2*nU+2*nl+njl) = df(:,2+2*nX+2*nU+nl+njl:1+2*nX+2*nU+2*nl+njl)*obj.options.lambda_mult;
      %       df(:,2+2*nX+2*nU+2*nl+2*njl:1+2*nX+2*nU+3*nl+2*njl) = df(:,2+2*nX+2*nU+2*nl+2*njl:1+2*nX+2*nU+3*nl+2*njl)*obj.options.Lambda_mult;
    end
    
    function [f,df] = dynamics_constraint_fun(obj,h,x0,x1,u0,u1,lambda0,lambda_jl0,lambda1,lambda_jl1,data0,data1)
      use_data = (nargin > 10);
      
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      nl = length(lambda0);
      njl = length(lambda_jl0);
      
      if use_data
        [xdot0,dxdot0] = constrained_dynamics(obj,x0(1:nq),x0(nq+1:end),u0,lambda0,lambda_jl0,data0);
        [xdot1,dxdot1] = constrained_dynamics(obj,x1(1:nq),x1(nq+1:end),u1,lambda1,lambda_jl1,data1);
      else
        [xdot0,dxdot0] = constrained_dynamics(obj,x0(1:nq),x0(nq+1:end),u0,lambda0,lambda_jl0);
        [xdot1,dxdot1] = constrained_dynamics(obj,x1(1:nq),x1(nq+1:end),u1,lambda1,lambda_jl1);
      end
      
      
      % cubic interpolation to get xcol and xdotcol, as well as
      % derivatives
      xcol = .5*(x0+x1) + h/8*(xdot0-xdot1);
      dxcol = [1/8*(xdot0-xdot1), (.5*eye(nX) + h/8*dxdot0(:,1:nX)), ...
        (.5*eye(nX) - h/8*dxdot1(:,1:nX)), h/8*dxdot0(:,nX+1:nX+nU) -h/8*dxdot1(:,nX+1:nX+nU), ...
        h/8*dxdot0(:,nX+nU+1:end), -h/8*dxdot1(:,nX+nU+1:end)];
      
      xdotcol = -1.5*(x0-x1)/h - .25*(xdot0+xdot1);
      dxdotcol = [1.5*(x0-x1)/h^2, (-1.5*eye(nX)/h - .25*dxdot0(:,1:nX)), ...
        (1.5*eye(nX)/h - .25*dxdot1(:,1:nX)), -.25*dxdot0(:,nX+1:nX+nU), -.25*dxdot1(:,nX+1:nX+nU),...
        -.25*dxdot0(:,nX+nU+1:end), -.25*dxdot1(:,nX+nU+1:end)];
      
      % evaluate xdot at xcol, using foh on control input
      [g,dgdxcol] = constrained_dynamics(obj,xcol(1:nq),xcol(nq+1:end),.5*(u0+u1),.5*(lambda0+lambda1),.5*(lambda_jl0+lambda_jl1));
      
      dg = dgdxcol(:,1:nX)*dxcol + [zeros(nX,1+2*nX), .5*dgdxcol(:,1+nX:nX+nU), .5*dgdxcol(:,1+nX:nX+nU),...
        .5*dgdxcol(:,1+nX+nU:end), .5*dgdxcol(:,1+nX+nU:end)];
      
      
      % constraint is the difference between the two
      f = xdotcol - g;
      df = dxdotcol - dg;
    end
    
    function [xdot,dxdot] = constrained_dynamics(obj,q,v,u,lambda,lambda_jl,data)
      use_data = (nargin > 6);
      
      %#ok<*MINV>
      nq = obj.plant.getNumPositions;
      nv = obj.plant.getNumVelocities;
      nu = obj.plant.getNumInputs;
      nl = length(lambda);
      njl = length(lambda_jl);
      
      assert(nq == nv,'Needs to be updated for nq != nv');
      
      [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(q,v);
      
      if use_data
        J = data.J;
        dJ = data.dJ;
      else
        [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
        J = zeros(nl,nq);
        dJ = zeros(nl*nq,nq);
        
        if obj.nC > 0
          J(1:1+obj.nD:end,:) = n;
          dJ(1:1+obj.nD:end,:) = dn;
          
          for j=1:length(D),
            J(1+j:1+obj.nD:end,:) = D{j};
            dJ(1+j:1+obj.nD:end,:) = dD{j};
          end
        end
      end
      [~,J_jl] = jointLimitConstraints(obj.plant,q);
      
      Hinv = inv(H);
      vdot = Hinv*(B*u - C + J'*lambda + J_jl'*lambda_jl);
      
      %dtau is [d/dq d/dv] of (B*u - C + J'*lambda + J_jl'*lambda_jl)
      dtau = -dC;
      
      if nu > 0        
        dtau = dtau + matGradMult(dB,u);
      end
      
      if obj.nC > 0
        dtau = dtau + [matGradMult(dJ,lambda,true), zeros(nv,nv)];
      end
      
      dv = [zeros(nq) eye(nv) zeros(nv,nu + nl + njl)];
      
      dvdot = [-Hinv*matGradMult(dH(:,1:nq),vdot) + Hinv*dtau(:,1:nq),...
        +Hinv*dtau(:,1+nq:end), Hinv*B, Hinv*J', Hinv*J_jl'];
      
      xdot = [v;vdot];
      dxdot = [dv;dvdot];
    end
    
    function [xtraj,utraj,ltraj,ljltraj,z,F,info] = solveTrajFromZ0(obj,z0)
      [z,F,info,infeasible_constraint_name] = obj.solve(z0);
      xtraj = reconstructStateTrajectory(obj,z);
      if nargout>1, utraj = reconstructInputTrajectory(obj,z); end
      
      t = [0; cumsum(z(obj.h_inds))];
      if obj.nC>0
        ltraj = PPTrajectory(foh(t,reshape(z([obj.lm_inds;obj.lp_inds;obj.L_inds]),[],obj.N)));
      else
        ltraj = [];
      end
      if obj.nJL>0
        ljltraj = PPTrajectory(foh(t,reshape(z(obj.ljl_inds),[],obj.N)));
      else
        ljltraj = [];
      end
    end
    
    function [xtraj,utraj,ltraj,ljltraj,z,F,info] = solveTraj(obj,t_init,traj_init)
      [xtraj,utraj,z,F,info] = solveTraj@DirectTrajectoryOptimization(obj,t_init,traj_init);
      t = [0; cumsum(z(obj.h_inds))];
      if obj.nC>0
        ltraj = PPTrajectory(foh(t,reshape(z([obj.lm_inds;obj.lp_inds;obj.L_inds]),[],obj.N)));
      else
        ltraj = [];
      end
      if obj.nJL>0
        ljltraj = PPTrajectory(foh(t,reshape(z(obj.ljl_inds),[],obj.N)));
      else
        ljltraj = [];
      end
    end
    
    function obj = setupVariables(obj,N)
      obj = setupVariables@DirectTrajectoryOptimization(obj,N);
      obj.nC = obj.plant.getNumContactPairs;
      [~,normal,d] = obj.plant.contactConstraints(zeros(obj.plant.getNumPositions,1));
      obj.nD = 2*length(d);
      assert(size(normal,2) == obj.nC); % just a double check
      
      nLambda_inds = obj.nC*(obj.nD+1);
      nContactForces = nLambda_inds*3 + 2*obj.nC;
      
      obj.lp_inds = reshape(obj.num_vars + (1:N*nLambda_inds),nLambda_inds,N);
      obj = obj.addDecisionVariable(N * nLambda_inds);
      
      obj.lm_inds = reshape(obj.num_vars + (1:N*nLambda_inds),nLambda_inds,N);
      obj = obj.addDecisionVariable(N * nLambda_inds);
      
      obj.L_inds = reshape(obj.num_vars + (1:N*nLambda_inds),nLambda_inds,N);
      obj = obj.addDecisionVariable(N * nLambda_inds);
      
      obj.gammap_inds = reshape(obj.num_vars + (1:N*obj.nC),obj.nC,N);
      obj = obj.addDecisionVariable(N * obj.nC);
      
      obj.gammam_inds = reshape(obj.num_vars + (1:N*obj.nC),obj.nC,N);
      obj = obj.addDecisionVariable(N * obj.nC);
      
      obj.dgamma_inds = reshape(obj.num_vars + (1:N*obj.nC),obj.nC,N);
      obj = obj.addDecisionVariable(N * obj.nC);
      
      %       obj.Gamma_inds = reshape(obj.num_vars + (1:N*obj.nC),obj.nC,N);
      %       obj = obj.addDecisionVariable(N * obj.nC);
      
      obj.phi_inds = reshape(obj.num_vars + (1:N*obj.nC),obj.nC,N);
      obj = obj.addDecisionVariable(N * obj.nC);
      
      obj.eta_inds = reshape(obj.num_vars + (1:N*obj.nC),obj.nC,N);
      obj = obj.addDecisionVariable(N * obj.nC);
      
      obj.nJL = obj.plant.getNumJointLimitConstraints();
      assert(obj.nJL == 0,'Not implemented yet--missing the impulses at least');
      
      obj.ljlp_inds = reshape(obj.num_vars + (1:N * obj.nJL),obj.nJL,N);
      obj = obj.addDecisionVariable(N * obj.nJL);
      
      obj.ljlm_inds = reshape(obj.num_vars + (1:N * obj.nJL),obj.nJL,N);
      obj = obj.addDecisionVariable(N * obj.nJL);
      
      obj.Ljl_inds = reshape(obj.num_vars + (1:N * obj.nJL),obj.nJL,N);
      obj = obj.addDecisionVariable(N * obj.nJL);
      
      obj.etajl_inds = reshape(obj.num_vars + (1:N * obj.nJL),obj.nJL,N);
      obj = obj.addDecisionVariable(N * obj.nJL);
      
      % joint limit constraints
      [jl_lb,jl_ub] = obj.plant.getJointLimits();
      obj.jl_lb_ind = find(jl_lb ~= -inf);
      obj.jl_ub_ind = find(jl_ub ~= inf);
    end
    
    % evaluates the initial trajectories at the sampled times and
    % constructs the nominal z0. Overwrite to implement in a different
    % manner
    function z0 = getInitialVars(obj,t_init,traj_init)
      if isscalar(t_init)
        t_init = linspace(0,t_init,obj.N);
      elseif length(t_init) ~= obj.N
        error('The initial sample times must have the same length as property N')
      end
      z0 = zeros(obj.num_vars,1);
      z0(obj.h_inds) = diff(t_init);
      
      if isfield(traj_init,'u') && obj.plant.getNumInputs > 0
        z0(obj.u_inds) = traj_init.u.eval(t_init);
      else
        nU = getNumInputs(obj.plant);
        z0(obj.u_inds) = 0.01*randn(nU,obj.N);
      end
      
      if isfield(traj_init,'x')
        z0(obj.x_inds) = traj_init.x.eval(t_init);
      else
        if ~isfield(traj_init,'u')
          traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,reshape(z0(obj.u_inds),nU,obj.N))),getInputFrame(obj.plant));
        end
        
        % todo: if x0 and xf are equality constrained, then initialize with
        % a straight line from x0 to xf (this was the previous behavior)
        
        %simulate
        sys_ol = cascade(traj_init.u,obj.plant);
        [~,x_sim] = sys_ol.simulate([t_init(1) t_init(end)]);
        z0(obj.x_inds) = x_sim.eval(t_init);
      end
      
      if obj.nC > 0
        if isfield(traj_init,'lm')
          z0(obj.lm_inds) = traj_init.lm.eval(t_init);
        end
        if isfield(traj_init,'lp')
          z0(obj.lp_inds) = traj_init.lp.eval(t_init);
        end
        if isfield(traj_init,'L')
          z0(obj.L_inds) = traj_init.L.eval(t_init);
        end
      end
      if obj.nJL > 0
        if isfield(traj_init,'ljl')
          z0(obj.ljl_inds) = traj_init.ljl.eval(t_init);
        else
          %           z0(obj.ljl_inds) = 0;
        end
      end
    end
    
    function obj = addRunningCost(obj,running_cost_function)
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      running_cost = FunctionHandleObjective(1+nX+nU,running_cost_function);
      for i=1:obj.N-1,
        obj = obj.addCost(running_cost,{obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)});
      end
    end
    
  end
end