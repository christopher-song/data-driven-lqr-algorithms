%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for quadrature with chebyshev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [kron_x_diff,x_diff,x_x,x_u,x_int,u_int] = chebyshev(dim_x, dim_u, x, u, dt, wts)
    
    samples = size(x,2);
    kron_xx= zeros(dim_x*dim_x,samples);
    kron_xu= zeros(dim_x*dim_u,samples);
    for j = 1:samples % compute kronecker products
        kron_xx(:,j) = kron(x(:,j), x(:,j));
        kron_xu(:,j) = kron(x(:,j), u(:,j));
    end
    kron_x_diff = kron_xx(:,end) - kron_xx(:,1); % get outputs
    x_diff = x(:,end) - x(:,1);
    x_x = 0.5*dt*(wts*(kron_xx.')).';
    x_u = 0.5*dt*(wts*(kron_xu.')).';
    x_int = 0.5*dt*(wts*(x.')).';
    u_int = 0.5*dt*(wts*(u.')).';
end