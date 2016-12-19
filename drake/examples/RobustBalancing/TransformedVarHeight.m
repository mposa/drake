classdef TransformedVarHeight < NStepCapturabilitySOSSystem
  % point foot, variable height, constant angular momentum
  % control input is (state-dependent scaling away from) force exerted
  % along vector from foot to CoM
  %
  % y = [x;xdot;1/z;zdot/z]
  % inputs are u1=foot COP
  %            u2 = z-accel (in g's)
  %            linearized so that xddot is missing a (-1/z*u1*u2*g) term
  % 
  % see mposa notes page 137
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    foot_radius
    f_range
  end
  
  methods
    function obj = TransformedVarHeight(g, z_nom, step_max, step_time, foot_radius, f_range)
      obj@NStepCapturabilitySOSSystem(4, 2, 1);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.f_range = f_range;
      obj.foot_radius = foot_radius;
    end
    
    function [xdot,dxdot] = dynamics(obj, t, x, u)
      [f,g] = controlAffineDynamics(obj, t, x);
      
%       f(3:4) = f(3:4)*0;
%       g(3:4,:) = g(3:4,:)*0;
%       f = subs(f,x(3:4),[0;0]);
%       g = subs(g,x(3:4),[0;0]);
      xdot = f + g*u;
      if nargout > 1
        dxdot = double(subs(obj.dxdot,obj.mssvars,[t;x;u]));
      end
    end
    
    function [f,g] = controlAffineDynamics(obj, t, y)
      qx = (y(1)+y(2))*sqrt(obj.z_nom/obj.g)/2;
      vx = y(2) - sqrt(obj.g/obj.z_nom)*qx;
      
      z_inv = (1/obj.z_nom+y(3));
      zddot_f = 0;
      zddot_g = [0 obj.g];
      xddot_f = z_inv*(qx*obj.g);
      xddot_g = z_inv*[-obj.g qx*obj.g];
      
      f = [sqrt(obj.g/obj.z_nom)*vx - xddot_f;sqrt(obj.g/obj.z_nom)*vx + xddot_f;...
        -y(4)*z_inv;zddot_f*z_inv - z_inv*y(4)^2];
      g = [-xddot_g;xddot_g;0 0;zddot_g*z_inv]*diag([obj.foot_radius; obj.f_range]);
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes x position only
      xp = xm;
      xp(1:2) = xm(1:2) + sqrt(obj.g/obj.z_nom)*s*obj.step_max;
    end
    
    
    % f_z / (m * g) <= f_max
    % 
    % f / (m * g) = u * |q| <= f_max
    % (f / (m * g))^2 = u^2 * q' * q <= f_max^2
    % and u > 0
    function ret = inputLimits(obj, u, x)
      ret = [1-u(1)^2;1-u(2)^2];     
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
%       q = x(1:2);
%       z = q(2) + obj.z_nom;
      umin = -ones(2,1);
      umax = ones(2,1);
      A = [];
    end
    
    function [A,b,C,d] = unitBoxInputTransform(obj)
      A = eye(2);
      b = zeros(2,1);
      C = A;
      d = b;
    end
    
    function ret = resetInputLimits(obj, s)
      ret = obj.step_max^2 - s'*s;
    end
    
    function [smin,smax] = simpleResetInputLimits(obj,x)
      smin = -1;
      smax = 1;
    end

    function rp = stanceUpdate(obj,x,r,s)
      rp = r + [s;0]*obj.step_max;
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, ~)
      q = x(1 : 2);
      v = x(3 : 4);
      
      sub_vars = [x(3:end);t];
      sub_val = [0;0;0];
      plot_vars = [x(1);x(2)];
      
      r_ic = x(2);
%       r_ic = x(1) + x(2)*sqrt(obj.z_nom / obj.g);
      dN = lipmCaptureLimit(obj.T, obj.foot_radius, obj.step_max, obj.z_nom, obj.g, n); % theoretical max ICP distance
%       dN = .5;
      
      figure(1)
      contourSpotless([Wsol;h_X;r_ic'*r_ic],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[1 0 dN^2],{'b','r','g'});
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')
      
      figure(n*10+2)
      h=contourSpotless([Vsol;h_X;r_ic'*r_ic],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0 dN^2],{'k','r','g'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      set(h,'LineWidth',4)
      % 3d plot for t = 0, zdot = 0
%       hFig = figure(n * 10 + 3);
%       clf;
%       contourSpotless3D(subs(Vsol, [x(4); t], [0; 0]), [x(1); x(3); x(2)], 0, [R_diag(1); R_diag(3); R_diag(2)]);
%       xlabel('q_1'); ylabel('v_1'); zlabel('q_2');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo(hFig,[class(obj) '_V' num2str(n)]);
      end
    end
    
    function draw(obj,t,x)
      qx = (x(1)+x(2))*sqrt(obj.z_nom/obj.g)/2;
      qz = 1/(1/obj.z_nom+x(3)) - obj.z_nom;
      qtheta = x(4);
      if length(x) > 4
        x_stance = x(end-1);
      else
        x_stance = 0;
      end
      % draw line from origin to COM
      h=line(x_stance+[0;qx],[0;qz + obj.z_nom]);
      set(h,'LineWidth',3,'Color','red')
      h=line([-10 10],[0 0]);
      set(h,'LineWidth',5,'Color','black')
      radius = .1;
      rectangle('Position',[x_stance+qx-radius/2,qz+obj.z_nom-radius/2,radius,radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-.5 1.5])
      ylim([-.5 1.5])
      axis off
    end
  end
end
