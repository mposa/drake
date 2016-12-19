classdef ScaledFull2D < NStepCapturabilitySOSSystem
  % point foot, variable height and angular momentum
  % control input is ground reaction force
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    inertia_ratio; % ratio of inertia to mass I/m
    fx_range
    fz_range
    step_max; % max step distance
    foot_size
    T; % step time
  end
  
  methods
    function obj = ScaledFull2D(g, z_nom, step_max, step_time, inertia_ratio, fx_range, fz_range, foot_size)
      obj@NStepCapturabilitySOSSystem(6, 3, 1);
      obj.g = g;
      obj.inertia_ratio = inertia_ratio;
      obj.z_nom = z_nom;
      obj.fx_range = fx_range;
      obj.fz_range = fz_range;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.foot_size = foot_size;
    end
    
    % u_1 = fx/mg
    % u_2 = (fz-mg)/mg
    function [xdot,dxdot] = dynamics(obj, t, x, u)     
      q = x(1 : 3);
      q(2) = q(2) + obj.z_nom;
      v = x(4 : 6);

      vdot = [u(1)*obj.fx_range*obj.g; ...
              u(2)*obj.fz_range*obj.g; ...
              obj.g/obj.inertia_ratio*(-u(2)*q(1)*obj.fz_range + u(1)*q(2)*obj.fx_range - q(1) + u(3)*obj.foot_size)];
      xdot = [v; vdot];
      
      if nargout > 1
        dvdotdx = [zeros(2,6); obj.g/obj.inertia_ratio*[-u(2)*obj.fz_range-1, u(1)*obj.fx_range, 0] zeros(1,3)];
        dvdotdu = [obj.fx_range*obj.g 0 0; 0 obj.fz_range*obj.g 0; obj.g/obj.inertia_ratio* [q(2)*obj.fx_range -q(1)*obj.fz_range obj.foot_size]];
        dxdot = [zeros(3,1), zeros(3) eye(3) zeros(3); zeros(3,1) dvdotdx dvdotdu];
      end
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q(1) only
      % qp = qm - u
      qm = xm(1 : 3);
      vm = xm(4 : 6);
      xp = [qm - [1; 0; 0] * s; vm];
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
      A = eye(3);
      b = zeros(3,1);
      C = A;
      d = b;
    end       
    
    function ret = resetInputLimits(obj, s)
      ret = obj.step_max^2 - s'*s;
    end

    function [smin,smax] = simpleResetInputLimits(obj,x)
      smin = -obj.step_max;
      smax = obj.step_max;
    end

    function rp = stanceUpdate(obj,x,r,s)
      rp = r + [s;0];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, u)
      q = x(1:3);
      v = x(4:6);
      
      sub_vars = [q(2:3);v(2:3);t];
      sub_val = [0;0;0;0;0];
      plot_vars = [q(1);v(1)];
      
      if ~isempty(Wsol)
        figure(1)
        contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[1 0],{'b','r'});
        xlabel('q_1')
        ylabel('v_1')
        title('W(x)')
      end
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
      % 3d plot for t = 0, zdot = 0
      contour_inds = [1 4 2];
      sub_inds = setdiff(1:6,contour_inds);
      hFig = figure(n * 10 + 3);
      clf;
      contourSpotless3D(subs(Vsol,  [t;x(sub_inds)], [0;0;0;0]), x(contour_inds), 0, R_diag(contour_inds));
      xlabel('x_c_m'); ylabel('xdot_c_m'); zlabel('z_c_m');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo(hFig,[class(obj) '_V' num2str(n)]);
      end
    end
    
    function draw(obj,t,x)
      if length(x) == obj.num_states
        x_stance = 0;
      else
        x_stance = x(end-1);
      end
      
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
