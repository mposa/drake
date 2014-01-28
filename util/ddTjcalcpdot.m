function [ddTJdot] = ddTjcalcpdot(jcode,q,qd)
% time derivative of single joint kinematic transformation for planar models
% this is the kinematic analog to jcalcp from the spatial vector library


switch (jcode)
  case 1 % pin joint
    c = cos(q);
    s = sin(q);
    ddTJdot(:,1) = [s*qd; c*qd; 0; -c*qd; s*qd; 0; 0; 0; 0];
    ddTJdot(:,2) = [-c; s; 0; -s; -c; 0; 0; 0; 0];
    ddTJdot(:,3) = ddTJdot(:,2);
    ddTJdot(:,4) = zeros(9,1);
  case 2 % x-axis prismatic
    ddTJdot = zeros(9,4);
  case 3 % y-axis prismatic
    ddTJdot = zeros(9,4);
end
