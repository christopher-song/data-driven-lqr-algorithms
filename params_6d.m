function [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_6d()
    Q = zeros(6,6); 
    Q(1,1) = 100;
    Q(6,6) = 100; % state cost 
    R = eye(2); % input cost
    K0 = zeros(2,6); % initial control used for ADP iteration
    A = [ -0.4125,   -0.0248,    0.0741,    0.0089,         0,         0; 
         101.5873,   -7.2651,    2.7608,    2.8068,         0,         0;
           0.0704,    0.0085,   -0.0741,   -0.0089,         0,    0.0200;
           0.0878,    0.2672,         0,   -0.3674,    0.0044,    0.3962; 
          -1.8414,    0.0990,         0,         0,   -0.0343,   -0.0330; 
                0,         0,         0, -359.0000,  187.5364,  -87.0316];
    B = [-0.0042,    0.0064;
         -1.0360,    1.5849; 
          0.0042,         0;
          0.1261,         0;
               0,   -0.0168;
               0,         0];
    a = 1; % amplitude of noise
    f = [100*(rand(100,1) - 0.5*ones(100,1)), 100*(rand(100,1) - 0.5*ones(100,1))]; % frequencies of noise
    dim_x = size(A,1); % size of x
    dim_u = size(B,2); % size of u
    [K_true,P_true,poles] = lqr(A,B,Q,R); % compute the true K and P
end
