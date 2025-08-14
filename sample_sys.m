%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for sampled data by integrating system with ode45 on one interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_out, u_out, x_end] = sample_sys(t_range, dim_u, x0, A, B, K0, a, f)

    [t,x] = ode45(@sys, t_range, x0, odeset('AbsTol',1e-50));
    
    function dx = sys(t,x)
        e = [sum(a*sin(f(:,1)*t)); sum(a*sin(f(:,2)*t))]; % sinusoidal exploration noise
        e = e(1:dim_u,1); % crop e to the right size
        u = K0*x + e;
        dx = A*x + B*u;  
    end
    x_out = x.';
    e_out = [sum(a*sin(kron(t_range,f(:,1))),1); sum(a*sin(kron(t_range,f(:,2))),1)]; % sinusoidal exploration noise
    e_out = e_out(1:dim_u,:); % crop e to the right size
    u_out = (K0*x.') + e_out; 
    x_end = x(end,:).';
end

