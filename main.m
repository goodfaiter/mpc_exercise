clear all
close all
yalmip('clear')
clc

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
R = 0.0;
N = 100;
P = diag([2 2 2 1 0 0 0]);

% Controller Variable Initialization
X = sdpvar(nx,N+1); % state trajectory: x0,x1,...,xN (columns of X)
Uin = sdpvar(nu,N); % input trajectory: u0,...,u_{N-1} (columns of U)

% Initialize objective and constraints of the problem
cost = 0.0; const = [];

% Assemble MPC formulation
for i = 1:N
    
    ref = [1.0
    0.1745 * sin(sys.Ts * i)
    -0.1745
    pi/2
    0
    0
    0];   
    
    % cost
    if( i < N )
        cost = cost + (X(:,i+1)-ref)'*Q*(X(:,i+1)-ref) + (Uin(:,i)-pinv(B)*(ref-A*ref))'*R*(Uin(:,i)-pinv(B)*(ref-A*ref));
    else
        cost = cost + (X(:,N+1)-ref)'*P*(X(:,N+1)-ref) + (Uin(:,N)-pinv(B)*(ref-A*ref))'*R*(Uin(:,N)-pinv(B)*(ref-A*ref));
    end
    
    % model
    const = [const, (X(:,i+1)-ref) == A*(X(:,i)-ref) + B*(Uin(:,i)-pinv(B)*(ref-A*ref))];
    
    % bounds
    const = [const, Umin <= Uin(:,i) <= Umax];
    const = [const, Xmin <= X(:,i+1) <= Xmax];
end

x0 = [0
    0
    0
    0
    0
    0
    0];

% Solve and plot
options = sdpsettings('solver','quadprog');
innerController = optimizer(const, cost, options, X(:,1), Uin(:,1));
simQuad( sys_inner, innerController, 0, x0, 10.0);

%%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
fprintf('PART V - simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Slew Rate Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VI - Slew Rate Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Soft Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VII - Soft Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FORCES Pro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VIII - FORCES Pro...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
