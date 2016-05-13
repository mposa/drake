classdef SimplePolynomialSystem < PolynomialSystem
  %SIMPLYPOLYNOMIALSYSTEM Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    f
    x
  end
  
  methods
    function obj = SimplePolynomialSystem(num_states,f,x)
      obj = obj@PolynomialSystem(num_states,0,0,num_states,false,true,false);
      obj.f = f;
      obj.x = x;
    end
    
     function xdot = dynamicsRHS(obj,t,x,u)
       xdot = subs(obj.f,obj.x,x);
     end
  end
  
end

