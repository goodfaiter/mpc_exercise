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

runFirstPart = false;
runSecondPart = false;
runThirdPart = true;

%%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System Linearization
% Akon = Ac(6:end,6:end);
% Bkon = Bc(6:end,:);
%
% [nx, nu] = size(Bkon);
%
% system_c = ss(Akon,Bkon,[],[]);
% system_d = c2d(system_c,sys.Ts);
%
% A = system_d.A;
% B = system_d.B;

A = sys.A;
B = sys.B;
[nx, nu] = size(B);

% Innerloop System Model
sys_inner = LTISystem('A', A, 'B', B, 'Ts', sys.Ts);

%%
% TODO: Run this onces properly and save the feasible set!
% lti = LTISystem('A', A, 'B', B);
%
% lti.x.min = [-0.1
%     degtorad(-1)
%     degtorad(-1)
%     degtorad(-180)
%     degtorad(-15)
%     degtorad(-15)
%     degtorad(-60)];
%
% lti.x.max = [0.1
%     degtorad(1)
%     degtorad(1)
%     degtorad(180)
%     degtorad(15)
%     degtorad(15)
%     degtorad(60)];
%
% lti.u.min = [0.0
%             0.0
%             0.0
%             0.0] - us;
% lti.u.max = [1.0
%             1.0
%             1.0
%             1.0] - us;
%
% InvSet = lti.invariantSet();
%
% A_f_in = InvSet.A;
% b_f_in = InvSet.b;
% for i = 1:size(A_f_in,1)
%     index = find(A_f_in(i,:));
%     if(A_f_in(i,index)<0)
%         X_fmin = [X_fmin b_f_in(i)/A_f_in(i,index)];
%     else
%         X_fmax = [X_fmax b_f_in(i)/A_f_in(i,index)];
%     end
% end
% X_fmin = X_fmin';
% X_fmax = X_fmax';
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

% MPC data
Q = diag([2 30 2 1 0 0 0]);
R = 0.1;
N = 19;
P = diag([2 30 2 1 0 0 0]);

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
    simQuad( sys_inner, innerController, 0, x0, 2.0);
end

%%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC data
Q = diag([2 2 2 1 0 0 0]);
R = 0.01;
N = 100;
P = diag([2 2 2 1 0 0 0]);

if (runSecondPart == true)
    
    % Controller Variable Initialization
    X = sdpvar(nx,N+1); % state trajectory: x0,x1,...,xN (columns of X)
    Uin = sdpvar(nu,N); % input trajectory: u0,...,u_{N-1} (columns of U)
    Ref = sdpvar(4,1);
    

    % Initialize objective and constraints of the problem
    cost = 0.0; const = [];
    
    % Delta-Formulation for tracking
%     This calculation gives the steady state solution to be:
%     xs = ref, us = 0;
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
    simQuad( sys_inner, innerController, 0, x0, T, r);

end

%%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Did this with innerController from Part II.

%%%%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC data
Q = diag([5 20 20 1 0 0 0]);
R = 0.01;
N = 100;
P = diag([5 20 20 1 0 0 0]);
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
        
        % cost
        if( i < N )
            cost = cost + (X(:,i+1)-ref)'*Q*(X(:,i+1)-ref) + Uin(:,i)'*R*Uin(:,i);
        else
            cost = cost + (X(:,N+1)-ref)'*P*(X(:,N+1)-ref) + Uin(:,N)'*R*Uin(:,N);
        end
        
        % model
        const = [const, X(:,i+1)-ref == A*X(:,i)-ref + B*Uin(:,i) + d(:,i)];
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
    r = [0
        0
        0
        0];
    
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
    simQuad( sys_inner, innerController, 0, x0, T, r, filter);
end

%%%%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
fprintf('PART V - simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Slew Rate Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VI - Slew Rate Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC data
Q = diag([2 2 2 1 0 0 0]);
R = 0.01;
N = 100;
P = diag([2 2 2 1 0 0 0]);

if (runSecondPart == true)
    
    % Controller Variable Initialization
    X = sdpvar(nx,N+1); % state trajectory: x0,x1,...,xN (columns of X)
    Uin = sdpvar(nu,N); % input trajectory: u0,...,u_{N-1} (columns of U)
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
    simQuad( sys_inner, innerController, 0, x0, T, r);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Soft Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VII - Soft Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FORCES Pro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VIII - FORCES Pro...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
