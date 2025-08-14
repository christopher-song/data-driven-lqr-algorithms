function [Q, R, K0, A, B, a, f, dim_x, dim_u, K_true, P_true] = params_2d()
    Q = eye(2); % state cost
    R = 1; % input cost
    K0 = [0,0]; % initial control used for ADP iteration
    A = [-1,-1; 1,-2];
    B = [1;1];
    a = 1; % amplitude of noise
    f = [5, 0]; % frequency of noise
    dim_x = size(A,1); % size of x
    dim_u = size(B,2); % size of u
    [K_true,P_true,poles] = lqr(A,B,Q,R); % compute the true K and P
end
