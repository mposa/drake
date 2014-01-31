classdef BallController < DrakeSystem
  %BALLCONTROLLER Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    K
    x0
  end
  
  methods
    function obj = BallController(p,x0,K)
      obj = obj@DrakeSystem(0,0,18,6);
      obj = obj.setInputFrame(p.getStateFrame);
      obj = obj.setOutputFrame(p.getInputFrame);
      obj.x0 = x0;
      obj.K = K;
    end
    
    function u = output(obj,t,~,x)
      u = -obj.K*(x - obj.x0);
    end
  end
  
end

