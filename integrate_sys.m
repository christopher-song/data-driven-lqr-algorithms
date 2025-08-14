%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for interval data by integrating system with ode45 on one interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kron_x_diff, x_diff, x_x, x_u, x_int, u_int, x_end] = integrate_sys(t_range, dim_x, dim_u, x0, A, B, K0, a, f)

    [t,z] = ode45(@sys_wrapper, t_range, [x0; zeros(dim_x*(dim_x+dim_u) + dim_x + dim_u,1)], odeset('AbsTol',1e-50));
    
    function dz = sys_wrapper(t,z)
        e = [sum(a*sin(f(:,1)*t)); sum(a*sin(f(:,2)*t))]; % sinusoidal exploration noise
        e = e(1:dim_u,1); % crop e to the right size
        x = z(1:dim_x); % the state of the system
        u = K0*x + e;
        % dimension of dx must match dim(x)!
        dx = A*x + B*u;  
        xx = kron(x,x);
        xu = kron(x,u);
        dz = [dx; xx; xu; x; u];
    end

    kron_x_diff = kron(z(end,1:dim_x).', z(end,1:dim_x).') - kron(z(1,1:dim_x).', z(1,1:dim_x).');
    x_diff = z(end,1:dim_x).' - z(1,1:dim_x).';
    x_x = z(end, dim_x+1: dim_x + dim_x*dim_x).';
    x_u = z(end, dim_x + dim_x*dim_x + 1: dim_x + dim_x*dim_x +dim_x*dim_u).';
    x_int = z(end, dim_x + dim_x*dim_x +dim_x*dim_u + 1: dim_x + dim_x + dim_x*dim_x + dim_x*dim_u).';
    u_int = z(end, dim_x + dim_x + dim_x*dim_x + dim_x*dim_u + 1: end).';
    x_end = z(end, 1:dim_x).';

end

% here is how z is allocated: 
% parts of each row of z with corresponding lengths:
% [     x,          xx,             xu,          x,          u]
% [dim(x),    dim(x)^2,  dim(x)*dim(u),     dim(x),     dim(u)]