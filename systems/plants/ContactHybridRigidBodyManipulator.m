classdef ContactHybridRigidBodyManipulator < HybridDrakeSystem
  properties
    plant
  end
  methods
    function obj = ContactHybridRigidBodyManipulator(p,constraint_ind_sequence)
      obj = obj@HybridDrakeSystem(p.getNumInputs,p.getNumStates);
      obj = setInputFrame(obj,p.getInputFrame);
      obj = setOutputFrame(obj,p.getOutputFrame);
      obj.plant = p;
      
      constraint_ind_sequence = constraint_ind_sequence(:);
      
      M = length(constraint_ind_sequence);
      % double up to do n and j
      for i=1:M,
        expanded_inds = sort([2*constraint_ind_sequence;2*constraint_ind_sequence+1]);
        obj = obj.addMode(ConstrainedRigidBodyManipulator(p,expanded_inds));
      end
      
      for i=1:M-1,
        obj = obj.addTransition(i,@(t,x,u) obj.guard_fun(t,x,u,constraint_ind_sequence{i},constraint_ind_sequence{i+1}),...
          @(mode,t,xm,u) obj.transition_fun(mode,t,xm,u,constraint_ind_sequence{i},constraint_ind_sequence{i+1}),...
          false,true,i+1);
        
      end
    end
    function [g,dg] = guard_fun(obj,t,x,u,ind_in,ind_out)
      i_transition = ind_out(ind_out > ind_in);
      if any(i_transition)
        q = x(1:obj.plant.getNumPositions);
        [phi,normal,d,xA,xB,idxA,idxB,mu,n] = obj.plant.contactConstraints(q);
        assert(sum(i_transition)  == 1)
        g = phi(i_transition);
        dg = [0 n(i_transition,:) zeros(1,length(u))];
      else
        g = 0;
        dg = zeros(1,1+length(x)+length(u));
      end
    end
    
    function [xp,mode,status,dxp] = transition_fun(obj,mode,t,xm,u,ind_in,ind_out)
      i_transition = ind_out(ind_out > ind_in);
      nX = length(xm);
      nU = length(u);
      
      if any(i_transition)
        n_active = length(i_out);
        
        q = xm(1:obj.plant.getNumPositions);
        qd = xm(obj.plant.getNumPositions+1:end);
        
        [H,C,B] = obj.plant.manipulatorDynamics(q,qd);
        
        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q);
        
        J = zeros(n_active*length(D),size(n,2));
        J(1:n_active,:) = n(i_out,:);
        for i=1:length(D),
          J(n_active*i+1:n_active*(i+1),:) = D{i}(i_out,:);
        end        
        
        qdp = qd -inv(H)*J'*inv(J*inv(H)*J')*J*qd;
        
        dqdpdqdm = []; % ugh, fine.
      else
        xp = xm;
        dxp = [zeros(1,nX) eye(nX) zeros(nU,nX)];
      end
    end
  end
end