classdef BallPlant < PlanarRigidBodyManipulator
  %BALLPLANT Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    phi_mss
    n_mss
    d_mss
    q_mss
    
  end
  
  methods
    function obj = BallPlant
      obj = obj@PlanarRigidBodyManipulator('ball.urdf');
      obj.num_contacts = 3;
      
      R = [.1;.05;.05];
      obj.q_mss = msspoly('q',9);
      x = obj.q_mss(1:2:6);
      z = obj.q_mss(2:2:6);
      r_perp_1_2 = [x(2) - x(1); z(2) - z(1)];
      r_perp_2_3 = [x(3) - x(2); z(3) - z(2)];
      r_perp_3_1 = [x(1) - x(3); z(1) - z(3)];
      
      r_fric_1_2 = [-r_perp_1_2(2);r_perp_1_2(1)];
      r_fric_2_3 = [-r_perp_2_3(2);r_perp_2_3(1)];
      r_fric_3_1 = [-r_perp_3_1(2);r_perp_3_1(1)];
      
      obj.phi_mss = [r_perp_1_2'*r_perp_1_2 - (R(1) + R(2))^2;...
        r_perp_2_3'*r_perp_2_3 - (R(2) + R(3))^2;...
        r_perp_3_1'*r_perp_3_1 - (R(3) + R(1))^2];
      
      obj.n_mss = diff(obj.phi_mss,obj.q_mss);
      
      obj.d_mss = [[r_fric_1_2; R(1)*r_fric_1_2'*r_fric_1_2; -r_fric_1_2; R(2)*r_fric_1_2'*r_fric_1_2;zeros(3,1)],...
        [zeros(3,1); r_fric_2_3; R(2)*r_fric_2_3'*r_fric_2_3; -r_fric_2_3; R(3)*r_fric_2_3'*r_fric_2_3],...
        [-r_fric_3_1; R(1)*r_fric_3_1'*r_fric_3_1; zeros(3,1); r_fric_3_1; R(3)*r_fric_3_1'*r_fric_3_1]]';
      
      
      % Transform from x1,y1,th1,x2,y2,th2,... to
      % x1,y1,x2,y2,...th1,th2,...
      T = zeros(9);
      T(1:2,1:2) = eye(2);
      T(3:4,4:5) = eye(2);
      T(5:6,7:8) = eye(2);
      T(7:9,[3 6 9]) = eye(3);
      
%       obj.n_mss = obj.n_mss*T;
      obj.d_mss = obj.d_mss*T';
    end
    
    function [phi,n,D,mu] = contactConstraints(obj, q)
      mu = 1*ones(3,1);
%       phi = double(subs(obj.phi_mss,obj.q_mss,q));
%       n = double(subs(obj.n_mss,obj.q_mss,q));
%       D{1} = double(subs(obj.d_mss,obj.q_mss,q));
%       D{2} = -D{1};
q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);
q5 = q(5);
q6 = q(6);
q7 = q(7);
q8 = q(8);
q9 = q(9);

phi = [
[  (-0.0225)+q1^2+q2^2+q3^2+(-2)*q3*q1+q4^2+(-2)*q4*q2  ]
[    (-0.01)+q3^2+q4^2+q5^2+(-2)*q5*q3+q6^2+(-2)*q6*q4  ]
[  (-0.0225)+q1^2+q2^2+q5^2+(-2)*q5*q1+q6^2+(-2)*q6*q2  ]];

n = [[  (2)*q1+(-2)*q3  (2)*q2+(-2)*q4  (-2)*q1+(2)*q3  (-2)*q2+(2)*q4               0               0    0    0    0  ]
[               0               0  (2)*q3+(-2)*q5  (2)*q4+(-2)*q6  (-2)*q3+(2)*q5  (-2)*q4+(2)*q6    0    0    0  ]
[  (2)*q1+(-2)*q5  (2)*q2+(-2)*q6               0               0  (-2)*q1+(2)*q5  (-2)*q2+(2)*q6    0    0    0  ]];

D{1} = [[  q2-q4  -q1+q3  -q2+q4   q1-q3       0      0  (0.1)*q1^2+(0.1)*q2^2+(0.1)*q3^2+(-0.2)*q3*q1+(0.1)*q4^2+(-0.2)*q4*q2  (0.05)*q1^2+(0.05)*q2^2+(0.05)*q3^2+(-0.1)*q3*q1+(0.05)*q4^2+(-0.1)*q4*q2                                                                          0  ]
[      0       0   q4-q6  -q3+q5  -q4+q6  q3-q5                                                                      0  (0.05)*q3^2+(0.05)*q4^2+(0.05)*q5^2+(-0.1)*q5*q3+(0.05)*q6^2+(-0.1)*q6*q4  (0.05)*q3^2+(0.05)*q4^2+(0.05)*q5^2+(-0.1)*q5*q3+(0.05)*q6^2+(-0.1)*q6*q4  ]
[  q2-q6  -q1+q5       0       0  -q2+q6  q1-q5  (0.1)*q1^2+(0.1)*q2^2+(0.1)*q5^2+(-0.2)*q5*q1+(0.1)*q6^2+(-0.2)*q6*q2                                                                          0  (0.05)*q1^2+(0.05)*q2^2+(0.05)*q5^2+(-0.1)*q5*q1+(0.05)*q6^2+(-0.1)*q6*q2  ]];
D{2} = -D{1};
%       x1 = q(1:2);
%       x2 = q(3:4);
%       x3 = q(5:6);
%       
%       phi = [norm(x2 - x1) - R(1) - R(2);...
%         norm(x3 - x2) - R(3) - R(2);...
%         norm(x1 - x3) - R(1) - R(3)];
%       
%       n = zeros(3,9);
%       n(1,:) = [x1(1) - x2(1); x1(2) - x2(2); x2(1) - x1(1); x2(2) - x1(2);zeros(5,1)]/phi(1);
%       n(2,:) = [zeros(2,1);x2(1) - x3(1); x2(2) - x3(2);x3(1) - x2(1); x3(2) - x2(2);zeros(3,1)]/phi(2);
%       n(3,:) = [x1(1) - x3(1); x1(2) - x3(2);zeros(2,1);x3(1) - x1(1); x3(2) - x1(2);zeros(3,1)]/phi(3);

    end
    
    function obj = compile(obj)
      obj = compile@PlanarRigidBodyManipulator(obj);
      obj.num_contacts = 3;
    end
  end
  
end

