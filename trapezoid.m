%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for quadrature with trapezoids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [kron_x_diff,x_diff,x_x,x_u,x_int,u_int] = trapezoid(dim_x, dim_u, x, u, dt, nt)
    times = linspace(0,dt,nt);
    samples = size(x,2);
    kron_xx= zeros(dim_x*dim_x,samples);
    kron_xu= zeros(dim_x*dim_u,samples);
    for j = 1:samples % compute kronecker products
        kron_xx(:,j) = kron(x(:,j), x(:,j));
        kron_xu(:,j) = kron(x(:,j), u(:,j));
    end
    kron_x_diff = kron_xx(:,end) - kron_xx(:,1); % get outputs
    x_diff = x(:,end) - x(:,1);
    x_x = trapz(times,kron_xx,2);
    x_u = trapz(times,kron_xu,2);
    x_int = trapz(times,x,2);
    u_int = trapz(times,u,2);
end