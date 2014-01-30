%converting rimless_roa_smaller_full_features.m to spotless
megaclear
sos_option = 1;

switch sos_option
  case 1
    sos_fun = @spot_mosek_sos;
    withSOS_fun = @withSOS;
  case 2
    sos_fun = @spot_mosek_sdsos;
    withSOS_fun = @withSDSOS;
  case 3
    sos_fun = @spot_mosek_dsos;
%     sos_fun = @spot_gurobi_dsos;
    withSOS_fun = @withDSOS;
end

for loop_count=2:2,
%   save tmp loop_count
%   clear all
%   load tmp
  p = PlanarRigidBodyManipulator('SimpleRimlessWheel_varchange.urdf');
  nC = p.getNumContacts;
      
  %**************
  % list of options
  %**************
  
  iter = rem(loop_count,2);
  inner_loop_max = 1;%iter*5 + 1;
  
  for inner_loop=1:inner_loop_max,
%     prog = mssprog();
    prog = spotsosprog();
    V_degree = 4;
    convexV = 0;
    
    ep = 1e-3;
    
    useVmult = 0;  % has no effect if searching for V
    
    % [z;s;1-c;qd]
%     Ai = diag([203;1;1;1;1;1]);
      [prog, Ai_diag] = prog.newFree(6,1);
      Ai = diag(Ai_diag);
      prog = prog.withPos(Ai_diag);
      prog = prog.withEqs(sum(Ai_diag) - 208);
%     [prog,A_inner] = new(prog,6,'psd');
%     prog.eq = trace(A_inner) - 208;
    
    %**************
    
    g = 9.81;
    
    q = msspoly('q',4);
    x = q(1);
    z = q(2);
    s = q(3);
    c = q(4);
    
    qd = msspoly('v',3);
    xd = qd(1);
    zd = qd(2);
    pitchd = qd(3);
    
    
    lx = msspoly('l',nC);
    
    lzsq = ones(nC,1);
    
    
    prog = prog.withIndeterminate(q(2:end));
    prog = prog.withIndeterminate(qd);
    prog = prog.withIndeterminate(lx);
    
    sos_eqns = [];
    
    % params = [];
    freevars = {};
    allfreevars = [x;z;xd;zd;s;c;pitchd;lx];
    
    switch iter
      case 0
        load iter_1_data
        
        [prog,d] = prog.newFree(1);
        
%         [prog,Ao] = prog.newPSD(6);
%         [prog,~] = prog.withEqs(trace(Ao) - 8.1);
        Ao = diag([2;2;3;.5;.5;.1]);
        
        %Create V
        [prog,V,coefv]= prog.newFreePoly(monomials([z;s;c;xd;zd;pitchd],0:V_degree));
        [prog,~] = prog.withEqs(subs(V,[z;s;c;xd;zd;pitchd],[0;0;1;0;0;0]));
        
      case 1
        if loop_count==1
          load(datapath('iter_0_data_default_spotless'))
          coefV = coefV*7/2;
          d = .2;
          % [z;s;1-c;qd]
          Ao = diag([2;2;3;.5;.5;.1]);
          V = reconstructPoly_spotless([z;s;c;xd;zd;pitchd],2,coefV);
        else
          load iter_0_data
          % do some semi-arbitrary scaling back of V
%           V = reconstructPoly_spotless([z;s;c;xd;zd;pitchd],V_degree,coefV)/.95;
        end
    end
    
    cpi8 = cos(pi/8);
    spi8 = sin(pi/8);
    rt2 = sqrt(2);
    
    
    f_impact=[0, 0;
      0, 0;
      0, 0;
      0, 0;
      lx(1), lx(2);
      lzsq(1)^2, lzsq(2)^2;
      4*lzsq(1)^2*(spi8*c + s*cpi8)- 4*(-lx(1))*(c*cpi8 - s*spi8), 4*lzsq(2)^2*(-spi8*c + s*cpi8)- 4*(-lx(2))*(c*cpi8 + s*spi8)];
    
    phi{1} = z + s*spi8 - c*cpi8 + cpi8;
    phi{2} = z - s*spi8 - c*cpi8 + cpi8;
    
    phidot{1} = zd + c*pitchd*spi8 + s*pitchd*cpi8;
    phidot{2} = zd - c*pitchd*spi8 + s*pitchd*cpi8;
    
    psi{1} = xd + c*pitchd*cpi8 - s*pitchd*spi8;
    psi{2} = xd + c*pitchd*cpi8 + s*pitchd*spi8;
    
    f_free=[xd;
      zd;
      c*pitchd;
      -s*pitchd;
      0;
      -g
      0];
    
    % E = (g*z + .5*zd^2) + .5*xd^2 + .5*.25*pitchd^2;
    
    Vdot_free = diff(V,[x;z;s;c;xd;zd;pitchd])*f_free;
    Vdot_impact_1 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,1);
    Vdot_impact_2 = diff(V,[x;z;s;c;xd;zd;pitchd])*f_impact(:,2);
    
    %% SOS functions
    x_vars = [q(2:end);qd];
    % Ball constraints
    [prog,rho] = prog.newFree(1);
    h_Bo = d - [z;s;1-c;qd]'*Ao*[z;s;1-c;qd];
    h_Bi = rho - [z;s;1-c;qd]'*Ai*[z;s;1-c;qd];
    
    % (1) -Vdot_free(x) >= 0 for x admissable in B_o
    % (2) -Vdot_impact_1(x,l) >= 0 for (x,l) admissable and x in B_o
    % (3) -Vdot_impact_2(x,l) >= 0 for (x,l) admissable and x in B_o
    % (4) V(x) >= 0 for x admissable and in B_o
    % (5) V(x) - 1 >= 0 for x on bdry(B_o)
    % (6) 1 - V(x) >= 0 for x in B_i
    
    sos_1 = -Vdot_free;
    sos_2 = -Vdot_impact_1;
    sos_3 = -Vdot_impact_2;
    sos_4 = V;
    sos_5 = (V - 1)*(1+x_vars'*x_vars);
%     sos_6 = -h_Bi*(1+x_vars'*x_vars);    
    sos_6 = 1-V;
    %% Add in constraints
    prog_bkp = prog;
    const_deg = V_degree;
    sig = {};
    coefsig = {};

    switch iter
      case 0
        
        % non-penetration (admissability of x)
        [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{1}, [x_vars], const_deg, sos_option);
        [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{2}, [x_vars], const_deg, sos_option);
        [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{1}*phi{2}, [x_vars], const_deg, sos_option);
        
        % trigonometric constraint
%         [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_1, 1 - s^2 - c^2, [x_vars], const_deg);
        
        sos_1 = sos_1 - h_Bo*reconstructPoly_spotless(x_vars,iter_1_mult.deg, iter_1_mult.coef);
        
        prog = withSOS_fun(prog,sos_1);
        
        %         Vsos_eqn = Vsos_eqn + (h-d)*reconstructPoly_spotless(freevars{1},Vlvlmult.degree,Vlvlmult.coef);
        
        % non-penetration (admissability of x)
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi{2}, [x_vars;lx(1)], const_deg, sos_option);
        
        % trigonometric constraint
%         [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, 1 - s^2 - c^2, [x_vars;lx(1)], const_deg);
        
        % Contact constraints (admissability of lambda)
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot{1}, [x_vars;lx(1)], const_deg, sos_option);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi{1}, [x_vars;lx(1)], const_deg, sos_option);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [x_vars;lx(1)], const_deg, sos_option);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi{1}, [x_vars;lx(1)], const_deg);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi{1}, [x_vars;lx(1)], const_deg);  %should this be psi^2?
               
        sos_2 = sos_2 - h_Bo*reconstructPoly_spotless([x_vars;lx(1)],iter_2_mult.deg, iter_2_mult.coef);
        
        prog = withSOS_fun(prog,sos_2);
        
        % non-penetration (admissability of x)
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi{1}, [x_vars;lx(2)], const_deg, sos_option);
        
        % trigonometric constraint
%         [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, 1 - s^2 - c^2, [x_vars;lx(2)], const_deg);
        
        % Contact constraints (admissability of lambda)
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot{2}, [x_vars;lx(2)], const_deg, sos_option);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi{2}, [x_vars;lx(2)], const_deg, sos_option);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg, sos_option);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi{2}, [x_vars;lx(2)], const_deg);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi{2}, [x_vars;lx(2)], const_deg);  %should this be psi^2?
        
        sos_3 = sos_3 - h_Bo*reconstructPoly_spotless([x_vars;lx(2)],iter_3_mult.deg, iter_3_mult.coef);
        
        prog = withSOS_fun(prog,sos_3);
        
        % non-penetration (admissability of x)
        [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{1}, x_vars, const_deg, sos_option);
        [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{2}, x_vars, const_deg, sos_option);
        [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{1}*phi{2}, x_vars, const_deg, sos_option);
        
        % trigonometric constraint
%         [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_4, 1 - s^2 - c^2, x_vars, const_deg);
        
        sos_4 = sos_4 - h_Bo*reconstructPoly_spotless(x_vars,iter_4_mult.deg, iter_4_mult.coef);
        prog = withSOS_fun(prog,sos_4);
        
        % non-penetration (admissability of x)
        [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{1}, x_vars, const_deg, sos_option);
        [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{2}, x_vars, const_deg, sos_option);
        [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{1}*phi{2}, x_vars, const_deg, sos_option);
        
        % trigonometric constraint
        [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, 1 - s^2 - c^2, x_vars, const_deg);
        
        sos_5 = sos_5 - h_Bo*reconstructPoly_spotless(x_vars,iter_5_mult.deg, iter_5_mult.coef);
        
        prog = withSOS_fun(prog,sos_5);
        
        % non-penetration (admissability of x)
        [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi{1}, x_vars, const_deg, sos_option);
        [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi{2}, x_vars, const_deg, sos_option);
        [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi{1}*phi{2}, x_vars, const_deg, sos_option);
        
        % trigonometric constraint
%         [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_6, 1 - s^2 - c^2, x_vars, const_deg);
        
        % Ball constraints
        sos_6 = sos_6 - (V-1)*reconstructPoly_spotless(x_vars,iter_6_mult_1.deg, iter_6_mult_1.coef);
        sos_6 = sos_6 - (V-1)*phi{1}*reconstructPoly_spotless(x_vars,iter_6_mult_2.deg, iter_6_mult_2.coef);
        sos_6 = sos_6 - (V-1)*phi{2}*reconstructPoly_spotless(x_vars,iter_6_mult_3.deg, iter_6_mult_3.coef);
        
        prog = withSOS_fun(prog,sos_6);
        cost = -rho;
        
      case 1
        cost = -rho;
        % non-penetration (admissability of x)
        [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{1}, [x_vars], const_deg, sos_option);
        [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{2}, [x_vars], const_deg, sos_option);
%         [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, phi{1}*phi{2}, [x_vars], const_deg, sos_option);
        
        % trigonometric constraint
        %             [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_1, 1 - s^2 - c^2, [x_vars], const_deg);
        
        [prog, sos_1, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_1, h_Bo, [x_vars], const_deg, sos_option);
        
        prog = withSOS_fun(prog,sos_1);
        
        % non-penetration (admissability of x)
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, phi{2}, [x_vars;lx(1)], const_deg, sos_option);
        
        % trigonometric constraint
        %             [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, 1 - s^2 - c^2, [x_vars;lx(1)], const_deg);
        
        % Contact constraints (admissability of lambda)
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -phidot{1}, [x_vars;lx(1)], const_deg, sos_option);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, -lx(1)*psi{1}, [x_vars;lx(1)], const_deg, sos_option);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, lzsq(1)^2 - lx(1)^2, [x_vars;lx(1)], const_deg, sos_option);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, phi{1}, [x_vars;lx(1)], const_deg);
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_2, (lzsq(1)^2 - lx(1)^2)*psi{1}, [x_vars;lx(1)], const_deg);  %should this be psi^2?
        
        [prog, sos_2, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_2, h_Bo, [x_vars;lx(1)], const_deg, sos_option);
        prog = withSOS_fun(prog,sos_2);
        
        % non-penetration (admissability of x)
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, phi{1}, [x_vars;lx(2)], const_deg, sos_option);
        
        % trigonometric constraint
        %             [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, 1 - s^2 - c^2, [x_vars;lx(2)], const_deg);
        
        % Contact constraints (admissability of lambda)
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -phidot{2}, [x_vars;lx(2)], const_deg, sos_option);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, -lx(2)*psi{2}, [x_vars;lx(2)], const_deg, sos_option);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, lzsq(2)^2 - lx(2)^2, [x_vars;lx(2)], const_deg, sos_option);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, phi{2}, [x_vars;lx(2)], const_deg);
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_3, (lzsq(2)^2 - lx(2)^2)*psi{2}, [x_vars;lx(2)], const_deg);  %should this be psi^2?
        
        [prog, sos_3, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_3, h_Bo, [x_vars;lx(2)], const_deg, sos_option);
        
        prog = withSOS_fun(prog,sos_3);
        
        % non-penetration (admissability of x)
        [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{1}, x_vars, const_deg, sos_option);
        [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{2}, x_vars, const_deg, sos_option);
%         [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, phi{1}*phi{2}, x_vars, const_deg, sos_option);
        
        % trigonometric constraint
        %             [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_4, 1 - s^2 - c^2, x_vars, const_deg);
        
        [prog, sos_4, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_4, h_Bo, x_vars, const_deg, sos_option);
        
        prog = withSOS_fun(prog,sos_4);
        
        % non-penetration (admissability of x)
        [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{1}, x_vars, const_deg, sos_option);
        [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{2}, x_vars, const_deg, sos_option);
%         [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_5, phi{1}*phi{2}, x_vars, const_deg, sos_option);
        
        % trigonometric constraint
        %             [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, 1 - s^2 - c^2, x_vars, const_deg);
        
        [prog, sos_5, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_5, h_Bo, x_vars, const_deg);
        
        prog = withSOS_fun(prog,sos_5);
        
        % non-penetration (admissability of x)
        [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi{1}, x_vars, const_deg, sos_option);
        [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi{2}, x_vars, const_deg, sos_option);
%         [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, phi{1}*phi{2}, x_vars, const_deg, sos_option);
        
        % trigonometric constraint
        %             [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_eq_sprocedure(prog, sos_6, 1 - s^2 - c^2, x_vars, const_deg);
        
        % Ball constraints
        [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, h_Bi, x_vars, const_deg, sos_option);
%         [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, V-1, x_vars, const_deg, sos_option);
%         [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, (V-1)*phi{1}, x_vars, const_deg, sos_option);
%         [prog, sos_6, sig{end+1}, coefsig{end+1}] = spotless_add_sprocedure(prog, sos_6, (V-1)*phi{2}, x_vars, const_deg, sos_option);
        prog = withSOS_fun(prog,sos_6);
    end
    
    options = spotprog.defaultOptions;
    options.verbose = true;
    options.scale_monomials = true;
    options.trig.enable = true;
    options.trig.sin = s;
    options.trig.cos = c;
    options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
    options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
%     options.solveroptions.method = 0;
%     options.solveroptions.presolve = 2;
    sol = prog.minimize(cost,sos_fun,options);

    
%% cleanup, reoptimize
backoff = true;
reoptimize = false;

if ~isequal(cost,0)
  if backoff
    costval = double(sol.eval(cost));
    prog = prog.withPos(costval + .01*abs(costval) - cost);
%     options = spotprog.defaultOptions;
%     options.scale_monomials = true;
%     options.verbose = true;
    sol = prog.minimize(0,sos_fun,options);
    display(sol.eval(cost))
  end
  
  if reoptimize
    for i=1:length(sol.gramMatrices),
      options.basis_scale{i} = eye(size(sol.gramMatrices{i},1));
    end
    
    for i=1:length(sol.gramMatrices),
      options.basis_scale{i} = full(reshape(clean(sol.eval(sol.gramMatrices{i})),size(sol.gramMatrices{i},1),[]))*options.basis_scale{i};
    end
    sol = prog.minimize(-rho,sos_fun,options);
  end
end
%% save results
    switch iter
      case 0
        d = double(sol.eval(d));
        coefV = double(sol.eval(coefv));
        vfn = sol.eval(V);
        Ao = double(sol.eval(Ao))
        Ai = double(sol.eval(Ai))
        R = double(sol.eval(rho))
        save iter_0_data_smaller d coefV Ao Ai R
        save(datapath(sprintf('roa_data%d',loop_count)),'d','coefV','R','Ao','Ai')
      case 1
        switch inner_loop
          case 1
            iter_1_mult.coef = double(sol.eval(coefsig{end}));
            iter_1_mult.deg = deg(sig{end},x_vars);
            save iter_1_data_1 iter_1_mult
          case 2
            iter_2_mult.coef = double(sol.eval(coefsig{end}));
            iter_2_mult.deg = deg(sig{end},x_vars);
            save iter_1_data_2 iter_2_mult
          case 3
            iter_3_mult.coef = double(sol.eval(coefsig{end}));
            iter_3_mult.deg = deg(sig{end},x_vars);
            save iter_1_data_3 iter_3_mult
          case 4
            iter_4_mult.coef = double(sol.eval(coefsig{end}));
            iter_4_mult.deg = deg(sig{end},x_vars);
            save iter_1_data_4 iter_4_mult
          case 5
            iter_5_mult.coef = double(sol.eval(coefsig{end}));
            iter_5_mult.deg = deg(sig{end},x_vars);
            save iter_1_data_5 iter_5_mult
          case 6
            iter_6_mult_1.coef = double(sol.eval(coefsig{end-2}));
            iter_6_mult_1.deg = deg(sig{end},x_vars);
            iter_6_mult_2.coef = double(sol.eval(coefsig{end-1}));
            iter_6_mult_2.deg = deg(sig{end},x_vars);
            iter_6_mult_3.coef = double(sol.eval(coefsig{end}));
            iter_6_mult_3.deg = deg(sig{end},x_vars);
            save iter_1_data_6 iter_6_mult_1 iter_6_mult_2 iter_6_mult_3
        end
    end
  end
end
