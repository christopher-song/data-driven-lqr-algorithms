function [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_3d()
    Q = eye(3); % state cost 
    R = 1; % input cost
    K0 = zeros(1,3); % initial control used for ADP iteration
    A = [   0,       1,       0;
            0,       0,       1;
         -0.1,    -0.5,    -0.7];
    B = [0; 0; 1];
    a = 2; % amplitude of noise
    f = [100,0]; % frequencies of noise
    dim_x = size(A,1); % size of x
    dim_u = size(B,2); % size of u
    [K_true,P_true,poles] = lqr(A,B,Q,R); % compute the true K and P
end
