classdef SwitchingController < DrakeSystem
  properties
    coeffs_B
    pows_B
    sos_plant
    nBaseStates
    nU
  end
  methods
    function obj = SwitchingController(plant,t_mss,x_mss,B,base_states)
      if nargin < 5
        base_states = true;
      end
      
      if base_states
        nBaseStates = 3;
      else
        nBaseStates = 0;
      end
      
      sos_plant = plant.sos_plant;
      obj = obj@DrakeSystem(0,0,sos_plant.num_states+nBaseStates,sos_plant.num_inputs);
      
      obj.nBaseStates = nBaseStates;
      obj.nU = length(B);
      for i=1:obj.nU,
        [pows_B,coeffs_B] = decomp_ordered(B(i),[t_mss;x_mss]);
        
        obj.coeffs_B{i} = coeffs_B;
        obj.pows_B{i} = pows_B;
      end
      

      obj = obj.setOutputFrame(plant.getInputFrame);
      obj = obj.setInputFrame(plant.getOutputFrame);
      obj.sos_plant = sos_plant;
    end
    
    function u = output(obj,t,~,x)
      if obj.nBaseStates > 0
        t = t - x(end-2);
        x = x(1:end-3);
      end
      u = zeros(obj.nU,1);
      [umin, umax, A] = obj.sos_plant.simpleInputLimits(x);
      for i=1:obj.nU,
        B = obj.coeffs_B{i}*prod(repmat([t;x]',length(obj.coeffs_B{i}),1).^obj.pows_B{i},2);
        if B < 0
          u(i) = umin(i);
        else
          u(i) = umax(i);
        end
      end
    end
  end
end