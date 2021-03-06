clear all
close all
yalmip('clear')
clc

%MPC Project

% adding the subfolders to the path
addpath(genpath('YALMIP-master'))
addpath(genpath('functions'))
addpath(genpath('data'))
addpath(genpath('FORCES_PRO_v1_6'))

% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outerController optimizer instance for the outer controller
load('quadData.mat')
outerController = getOuterController(Ac);
disp('Data successfully loaded')

%%%%%%%%%%%%%%%% ADD YOUR CODE BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

part = 5;
runFirstPart = false;
runSecondPart = false;  
runThirdPart = false;
runFourthPart = false;
runFifthPart = false;
invariantSet = false;

switch part
    case 0.1
        invariantSet = true;
    case 1
        runFirstPart = true;
    case 2
        runSecondPart = true;
    case 3
        runThirdPart = true;
    case 4
        runFourthPart = true;
    case 5
        runFifthPart = true;
    otherwise
        runFirstPart = false;
        runSecondPart = false;
        runThirdPart = false;
        runFourthPart = false;
        invariantSet = false;
        disp('Enter valid number')
end

%%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = sys.A;
B = sys.B;
[nx, nu] = size(B);

% Innerloop System Model
sys_inner = LTISystem('A', A, 'B', B, 'Ts', sys.Ts);

%%

% MPC data
Q = diag([5 100 100 1 0 0 0]);
R = 0.05*eye(4);
N = 20;
P = 100*diag([5 20 20 1 0 0 0]);

% TODO: Run this onces properly and save the feasible set!
if invariantSet == true
    
    % Invariant Terminal Set Computation
      
%     lti = LTISystem('A', A, 'B', B, 'Ts', sys.Ts);
%     
%     lti.x.min = [-0.1
%         degtorad(-1)
%         degtorad(-1)
%         degtorad(-180)
%         degtorad(-15)
%         degtorad(-15)
%         degtorad(-60)];
%     
%     lti.x.max = [0.1
%         degtorad(1)
%         degtorad(1)
%         degtorad(180)
%         degtorad(15)
%         degtorad(15)
%         degtorad(60)];
%     
%     lti.u.min = [0.0
%         0.0
%         0.0
%         0.0] - us;
%     lti.u.max = [1.0
%         1.0
%         1.0
%         1.0] - us;
%     
%     InvSet = lti.invariantSet();

    % Terminal Maximum Control Invariant Set Computation
    [P_inf,~,~] = dare(A,B,Q,R);
    F_inf = -inv(B'*P_inf*B + R)*B'*P_inf*A;
    
    lti = LTISystem('A', A+B*F_inf);

    A_poly = [-eye(7);-F_inf;eye(7);F_inf];      

    b_poly = [-0.1;
        -degtorad(-1);
        -degtorad(-1);
        -degtorad(-180);
        -degtorad(-15);
        -degtorad(-15);
        -degtorad(-60);
        -([0;0;0;0]-us);
        0.1;
        degtorad(1);
        degtorad(1);
        degtorad(180);
        degtorad(15);
        degtorad(15);
        degtorad(60);
        [1;1;1;1]-us];

    poly_set = Polyhedron( 'A', A_poly, 'b', b_poly);

    InvSet = lti.invariantSet('X', poly_set);
    % Result: Empty polyhedron in R^7
end
%%

% Constraint Initialization
Xmin = [-1.0
    degtorad(-10)
    degtorad(-10)
    degtorad(-180)
    degtorad(-15)
    degtorad(-15)
    degtorad(-60)];

Xmax = [1.0
    degtorad(10)
    degtorad(10)
    degtorad(180)
    degtorad(15)
    degtorad(15)
    degtorad(60)];

X_fmin = [-0.1
    degtorad(-1)
    degtorad(-1)
    degtorad(-180)
    degtorad(-15)
    degtorad(-15)
    degtorad(-60)];

X_fmax = [0.1
    degtorad(1)
    degtorad(1)
    degtorad(180)
    degtorad(15)
    degtorad(15)
    degtorad(60)];

Umin = [0.0
    0.0
    0.0
    0.0] - us;

Umax = [1.0
    1.0
    1.0
    1.0] - us;


if (runFirstPart == true)
    
    % Controller Variable Initialization
    X = sdpvar(nx,N+1); % state trajectory: x0,x1,...,xN (columns of X)
    Uin = sdpvar(nu,N); % input trajectory: u0,...,u_{N-1} (columns of U)
    
    % Initialize objective and constraints of the problem
    cost = 0.0; const = [];
    
    % Assemble MPC formulation
    for i = 1:N
        % cost
        if( i < N )
            cost = cost + X(:,i+1)'*Q*X(:,i+1) + Uin(:,i)'*R*Uin(:,i);
        else
            cost = cost + X(:,N+1)'*P*X(:,N+1) + Uin(:,N)'*R*Uin(:,N);
        end
        
        % model
        const = [const, X(:,i+1) == A*X(:,i) + B*Uin(:,i)];
        
        % bounds
        if( i < N )
            const = [const, Umin <= Uin(:,i) <= Umax];
            const = [const, Xmin <= X(:,i+1) <= Xmax];
        else
            const = [const, Umin <= Uin(:,N) <= Umax];
            const = [const, X_fmin <= X(:,N+1) <= X_fmax];
        end
    end
    
    T = 2;
    
    % Initial state
    x0 = [-1.0
        0.1745
        -0.1745
        0.8727
        0.0
        0.0
        0.0];
    
    % Solve and plot
    options = sdpsettings('solver','quadprog');
    innerController = optimizer(const, cost, options, X(:,1), Uin(:,1));
    simQuad( sys_inner, innerController, 0, x0, T);
end

%%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;

if (runSecondPart == true)
    
    % Controller Variable Initialization
    X = sdpvar(nx,N+1); % state trajectory: x0,x1,...,xN (columns of X)
    Uin = sdpvar(nu,N); % input trajectory: u0,...,u_{N-1} (columns of U)
    Ref = sdpvar(4,1);
    
    
    % Initialize objective and constraints of the problem
    cost = 0.0; const = [];
    
    % Delta-Formulation for tracking
    %     %This calculation gives the steady state solution to be:
    %     %xs = ref, us = 0;
    %     syms r1 r2 r3 r4;
    %     rr = [r1;r2;r3;r4];
    %     C = [eye(4) zeros(4,3)];
    %     temp = [eye(7)-A -B; C zeros(4)];
    %     temp2 = [zeros(7,1);rr];
    %     temp = [temp temp2]
    
    % Assemble MPC formulation
    for i = 1:N
        
        
        ref = [Ref(:,1)
            0
            0
            0];
        
        % Delta-Formulation for tracking
        X_delta_k = (X(:,i)-ref);
        X_delta_k_1 = (X(:,i+1)-ref);
        X_delta_N = (X(:,N+1)-ref);
        
        % cost
        if( i < N )
            cost = cost + X_delta_k_1'*Q*X_delta_k_1 + Uin(:,i)'*R*Uin(:,i);
        else
            cost = cost + X_delta_N'*P*X_delta_N + Uin(:,N)'*R*Uin(:,N);
        end
        
        % model
        const = [const, X_delta_k_1 == A*X_delta_k + B*Uin(:,i)];
        
        % bounds
        const = [const, Umin <= Uin(:,i) <= Umax];
        const = [const, Xmin-ref <= X_delta_k_1 <= Xmax-ref];
    end
    
    x0 = [0
        0
        0
        0
        0
        0
        0];
    
    % Part 5 Reference
    T = 10;
    r = [1.0
        0.1745
        -0.1745
        1.7453];
    
    % Part 6 Reference
    T = 10;
    steps = floor(T/sys.Ts);
    for step = 1:steps
        r(1,step) = 1;
        r(2,step) = 0.1745*sin(step*sys.Ts);
        r(3,step) = -0.1745*sin(step*sys.Ts);
        r(4,step) = 1.7453;
    end
    
    % Solve and plot
    options = sdpsettings('solver','quadprog');
    innerController = optimizer(const, cost, options, [X(:,1)' Ref(:,1)']', Uin(:,1));
    
    % Forces Optimization
    %     codeoptions = getOptions('internal_simpleMPC_solver_1');
    %     innerController = optimizerFORCES(const, cost, codeoptions, [X(:,1)' Ref(:,1)']', Uin(:,1), {'xinit'}, {'u0'});
    [xt ut t rt deltat] = simQuad( sys_inner, innerController, 0, x0, T, r);
    
end

%%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC Data
Q = diag([5 100 100 1 20 20 20]);
Ld = 1.0*diag([1 1.5 1.5 1 0.01 0.01 0.01]);

if runThirdPart == true
    % Controller Variable Initialization
    X = sdpvar(nx,N+1); % state trajectory: x0,x1,...,xN (columns of X)
    Uin = sdpvar(nu,N); % input trajectory: u0,...,u_{N-1} (columns of U)
    d = sdpvar(nx,N+1);
    Ref = sdpvar(4,1);
    
    % Initialize objective and constraints of the problem
    cost = 0.0; const = [];
    
    % Assemble MPC formulation
    for i = 1:N
        
        ref = [Ref(:,1)
            0
            0
            0];
        % Delta-Formulation for tracking
        X_delta_k = (X(:,i)-ref);
        X_delta_k_1 = (X(:,i+1)-ref);
        X_delta_N = (X(:,N+1)-ref);
        
        % cost
        if( i < N )
            cost = cost + X_delta_k_1'*Q*X_delta_k_1 + Uin(:,i)'*R*Uin(:,i);
        else
            cost = cost + X_delta_N'*P*X_delta_N + Uin(:,N)'*R*Uin(:,N);
        end
        
        % model
        const = [const, X_delta_k_1 == A*X_delta_k + B*Uin(:,i) + d(:,i)];
        const = [const, d(:,i+1) == d(:,i)];
        
        % bounds
        
         const = [const, Umin <= Uin(:,i) <= Umax];
         const = [const, Xmin-ref <= X(1:7,i+1)-ref <= Xmax-ref];
    end
    
    A_aug = [A eye(nx); zeros(7) eye(nx)];
    B_aug = [B; zeros(7,4)];
    C_aug = [eye(nx) eye(nx)];
    
    L = [eye(nx); Ld];
    Af = A_aug - L*C_aug;
    Bf = [B_aug L];
    
    filter = struct('Af', Af, 'Bf', Bf);
    
    T = 15;
    
    r = [0.8
        0.12
        -0.12
        pi/2];
    
%     r = [0
%         0
%         0
%         0];
    
    steps = floor(T/sys.Ts);
    for step = 1:steps
        r(1,step) = 0.8;
        r(2,step) = 0.12*sin(step*sys.Ts);
        r(3,step) = -0.12*sin(step*sys.Ts);
        r(4,step) = pi/2;
    end
   
    x0 = zeros(7,1);
    % Solve and plot
    options = sdpsettings('solver','quadprog');
    innerController = optimizer(const, cost, options, [X(:,1)' Ref(:,1)' d(:,1)']', Uin(:,1));
    
    % Forces Optimization
    %     codeoptions = getOptions('internal_simpleMPC_solver_1');
    %     innerController = optimizerFORCES(const, cost, codeoptions, [X(:,1)' Ref(:,1)' d(:,1)']', Uin(:,1), {'xinit'}, {'u0'});
    [xt ut t rt deltat] = simQuad( sys_inner, innerController, 0, x0, T, r, filter);
end

%%%%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
fprintf('PART V - simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Slew Rate Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VI - Slew Rate Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC data
N = 100;
% Slew rate contraint
delta = 0.3*[1 1 1 1]';

if (runFourthPart == true)
    % Controller Variable Initialization
    X = sdpvar(nx,N+1);
    Uin = sdpvar(nu,N);
    d = sdpvar(nx,N+1);
    Ref = sdpvar(4,1);
    u_prev = sdpvar(nu,1);
    
    % Initialize objective and constraints of the problem
    cost = 0.0; const = [];
    
    % Assemble MPC formulation
    for i = 1:N
        
        ref = [Ref(:,1)
            0
            0
            0];
        
        % Delta-Formulation for tracking
        X_delta_k = (X(:,i)-ref);
        X_delta_k_1 = (X(:,i+1)-ref);
        X_delta_N = (X(:,N+1)-ref);
        
        % cost
        if( i < N )
            cost = cost + X_delta_k'*Q*X_delta_k + Uin(:,i)'*R*Uin(:,i);
        else
            cost = cost + X_delta_N'*P*X_delta_N + Uin(:,N)'*R*Uin(:,N);
        end
        
        % model
        const = [const, X_delta_k_1 == A*X_delta_k + B*Uin(:,i) + d(:,i)];
        const = [const, d(:,i+1) == d(:,i)];
        
        % bounds
        const = [const, Umin <= Uin(:,i) <= Umax];
        const = [const, Xmin-ref <= X(1:7,i+1)-ref <= Xmax-ref];
        
        % Slew Contraints
        if i == 1
            const = [const, Uin(:,i) - u_prev <= delta];
            const = [const, Uin(:,i) - u_prev >= -delta];
        elseif (i < N )
            const = [const, Uin(:,i+1) - Uin(:,i) <= delta];
            const = [const, Uin(:,i+1) - Uin(:,i) >= -delta];
        end
    end
    
    A_aug = [A eye(nx); zeros(7) eye(nx)];
    B_aug = [B; zeros(7,4)];
    C_aug = [eye(nx) eye(nx)];
    
    L = [eye(nx); Ld];
    Af = A_aug - L*C_aug;
    Bf = [B_aug L];
    
    filter = struct('Af', Af, 'Bf', Bf);
    
    T = 10;
    
    r = [0.8
        0.12
        -0.12
        pi/2];
    
    x0 = zeros(7,1);
    
    % Solve and plot
    options = sdpsettings('solver','quadprog');
    innerController = optimizer(const, cost, options, [X(:,1)' Ref(:,1)' u_prev' d(:,1)']', Uin(:,1:2));
    [xt ut t rt deltat] = simQuad( sys_inner, innerController, 0, x0, T, r, filter, [], 2);
    
    figure(11); clf; grid on; hold on;
    for i = 1:size(ut,2)-1
        delta_graph(1:4,i) = ut(1:4,i+1) - ut(1:4,i);
    end
    plot(t(1:end-1), delta_graph,'LineWidth',1.1);
    ylabel('Delta of Input'); xlabel('s');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Soft Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VII - Soft Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC Data
N = 100;

% Soft Constraints
s = 2*diag([1 1 1 1]);
v = 0.25*[1 1 1 1]';

if (runFifthPart == true)
    % Controller Variable Initialization
    X = sdpvar(nx,N+1);
    Uin = sdpvar(nu,N);
    d = sdpvar(nx,N+1);
    Ref = sdpvar(4,1);
    u_prev = sdpvar(nu,1);
    epsilon = sdpvar(nu,N); % slack variable
    
    % Initialize objective and constraints of the problem
    cost = 0.0; const = [];
    
    % Assemble MPC formulation
    for i = 1:N
        
        ref = [Ref(:,1)
            0
            0
            0];
        
        % Delta-Formulation for tracking
        X_delta_k = (X(:,i)-ref);
        X_delta_k_1 = (X(:,i+1)-ref);
        X_delta_N = (X(:,N+1)-ref);
        
        % cost
        if( i < N )
            cost = cost + X_delta_k'*Q*X_delta_k + Uin(:,i)'*R*Uin(:,i) ...
                + v'*epsilon(:,i) + epsilon(:,i)'*s*epsilon(:,i);
        else
            cost = cost + X_delta_N'*P*X_delta_N + Uin(:,N)'*R*Uin(:,N);
        end
        
        % model
        const = [const, X_delta_k_1 == A*X_delta_k + B*Uin(:,i) + d(:,i)];
        %const = [const, X_delta_k_1 == A*X_delta_k + B*Uin(:,i)];
        const = [const, d(:,i+1) == d(:,i)];
        
        % bounds
        const = [const, Umin <= Uin(:,i) <= Umax];
        const = [const, Xmin-ref <= X(1:7,i+1)-ref <= Xmax-ref];
        
        % Slew Contraints
        if i == 1
            const = [const, Uin(:,i) - u_prev <= delta + epsilon(:,i)];
            const = [const, Uin(:,i) - u_prev >= -(delta + epsilon(:,i))];
        elseif (i < N )
            const = [const, Uin(:,i+1) - Uin(:,i) <= delta + epsilon(:,i)];
            const = [const, Uin(:,i+1) - Uin(:,i) >= -(delta + epsilon(:,i))];
        end
        
        %Soft Constrants
        const = [const, 0 <= epsilon(:,i)];
    end
    
    A_aug = [A eye(nx); zeros(7) eye(nx)];
    B_aug = [B; zeros(7,4)];
    C_aug = [eye(nx) eye(nx)];
    
    L = [eye(nx); Ld];
    Af = A_aug - L*C_aug;
    Bf = [B_aug L];
    
    filter = struct('Af', Af, 'Bf', Bf);
    
    T = 10;
    
    r = [0.8
        0.12
        -0.12
        pi/2];
    
    x0 = zeros(7,1);
    
    % Solve and plot
    options = sdpsettings('solver','quadprog');
    innerController = optimizer(const, cost, options, [X(:,1)' Ref(:,1)' u_prev' d(:,1)']', Uin(:,1:2));
    [xt ut t rt deltat] = simQuad( sys_inner, innerController, 0, x0, T, r, filter, [], 2);
    
    figure(11); clf; grid on; hold on;
    for i = 1:size(ut,2)-1
        delta_graph(1:4,i) = ut(1:4,i+1) - ut(1:4,i);
        delta_graph_max(1:4,i) = delta';
        delta_graph_min(1:4,i) = -delta';
    end
    plot(t(1:end-1), delta_graph,'LineWidth',1.1);
    plot(t(1:end-1), delta_graph_max,'r--');
    plot(t(1:end-1), delta_graph_min,'r--');
    ylabel('Delta of Input'); xlabel('s');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FORCES Pro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VIII - FORCES Pro...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
