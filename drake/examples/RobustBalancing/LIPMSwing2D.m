classdef LIPMSwing2D < NStepCapturabilitySOSSystem
  % point foot, variable height and angular momentum
  % control input is ground reaction force
  
  properties
    gravity; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    foot_radius
    swing_speed
  end
  
  methods
    function obj = LIPMSwing2D(g, z_nom, step_max, step_time, foot_radius,swing_speed)
      obj@NStepCapturabilitySOSSystem(3, 2, 0);
      obj.gravity = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;      
      obj.foot_radius = foot_radius;
      obj.swing_speed = swing_speed;
    end
    
    % x = [q;v]
    % q = [x_com;z_com;theta]
    function xdot = dynamics(obj, t, x, u)
      [f,g] = controlAffineDynamics(obj, t, x);
      xdot = f + g*u;
    end
    
    function [f,g] = controlAffineDynamics(obj,t,x)
      q = x(1);
      v = x(2);
      
      z_inv = 1/obj.z_nom;
      
      f_x_ddot = obj.gravity*z_inv*q(1);
      g_x_ddot = -obj.gravity*z_inv;
      
      f = [v;f_x_ddot];
      g = [0; g_x_ddot];
      g = g*diag([obj.foot_radius]);
      
      % add swing leg
      f = [f;0];
      g = [[g;0] [zeros(2,1);obj.swing_speed]];
    end    
    
    function [xp,constraint] = reset(obj, t, xm, s)
      qm = xm(1);
      vm = xm(2);
      x_swing = xm(3);
      qp = qm(1) - x_swing;
      vp = vm;
      x_swingp = -x_swing;
      xp = [qp;vp;x_swingp];
      constraint = [];
    end
    
    function [xm,constraint] = inverse_reset(obj, t, xp, s)
      qp = xp(1);
      vp = xp(2);
      x_swingp = xp(3);
      
      qm = qp(1)-x_swingp;
      vm = vp;
      x_swingm = -x_swingp;
      xm = [qm;vm;x_swingm];
      
      constraint = [];
    end
    
    function ret = inputLimits(obj, u, x)
      ret = 1 - u.*u;
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      umin = -ones(3,1);
      umax = ones(3,1);
      A = [];
    end
    
    function [A,b,C,d] = unitBoxInputTransform(obj)
      A = 1;
      b = 0;
      C = A;
      d = b;
    end       
    
    function ret = resetInputLimits(obj, s)
      ret = [];
    end

    function [smin,smax] = simpleResetInputLimits(obj,x)
      smin = -inf
      smax = inf;
    end

    function rp = stanceUpdate(obj,x,r,s)
      rp = r + obj.step_max*[s;0];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, u)
      q = x(1);
      v = x(2);
      
      sub_vars = [x(3);t];
      sub_val = [0;0;];
      plot_vars = [q(1);v(1)];
      
      if ~isempty(Wsol)
        figure(1)
        contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[1 0],{'b','r'});
        xlabel('q_1')
        ylabel('v_1')
        title('W(x)')
      end
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
    end
    
    function draw(obj,t,x)
      x_stance = x(end-1);
      % draw line from origin to COM
      h=line(x_stance+[0;x(1)],[0;x(2) + obj.z_nom]);
      set(h,'LineWidth',3,'Color','red')
      axis_major = .3;
      axis_minor = .2;
      theta = linspace(0,2*pi,100);
      x_ellipse = axis_minor*cos(theta);
      z_ellipse = axis_major*sin(theta);
      x_body = x_stance + x(1) + x_ellipse*cos(x(3)) + z_ellipse*sin(x(3));
      z_body = obj.z_nom + x(2) - x_ellipse*sin(x(3)) + z_ellipse*cos(x(3));
      patch(x_body,z_body,'k')
%       rectangle('Position',[x_stance+x(1)-radius/2,x(2)+obj.z_nom-radius/2,radius,2*radius],'Curvature',[1,1], 'FaceColor','k')
%       xlim([-3 3])
%       ylim([-.1 1.5])
      
            h=line([-10 10],[0 0]);
      set(h,'LineWidth',5,'Color','black')
      radius = .1;
      rectangle('Position',[x_stance+x(1)-radius/2,x(2)+obj.z_nom-radius/2,radius,radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-.5 1.5])
      ylim([-.5 1.5])
      axis off
    end
  end
end
