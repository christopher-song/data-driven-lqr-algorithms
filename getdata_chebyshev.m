%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for getting chebyshev approximation from sampled data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kron_x_diffs, x_diffs, xxs, xus, xs, us] = getdata_chebyshev(dim_x, dim_u, intervals, dt, nt, x0, A, B, K0, a, f)
  
    % arrays for data from ADP
    kron_x_diffs = zeros(intervals,dim_x*dim_x);
    xxs = zeros(intervals,dim_x*dim_x);
    xus = zeros(intervals,dim_x*dim_u);
    
    % arrays for data from system ID
    x_diffs = zeros(intervals, dim_x);
    xs = zeros(intervals,dim_x);
    us = zeros(intervals,dim_u);

    [pts,wts] = chebpts(nt); % get chebyshev points and weights
    pts = ((0.5*(dt))*pts + (0.5*(dt))*ones(nt,1)).'; % scale points to our interval
   
    for i = 1:intervals
        %interval = i % print interval number
        [x, u, x_end] = sample_sys((i-1)*dt*ones(1,nt) + pts, dim_u, x0, A, B, K0, a, f); % integrate
        % do quadrature
        [kron_x_diff,x_diff,x_x,x_u,x_int,u_int] = chebyshev (dim_x, dim_u, x,u,dt,wts);
        % save data
        kron_x_diffs(i,:) = kron_x_diff.';
        x_diffs(i,:) = x_diff.';
        xxs(i,:) = x_x.';
        xus(i,:) = x_u.';
        xs(i,:) = x_int.';
        us(i,:) = u_int.';
        if dim_x == 2
            x0 = x_end; % update x0
        elseif dim_x == 3
            x0 = 10*(rand(3,1)-0.5*ones(3,1)); % update x0
        else % the dim is 6
            x0 = (rand(6,1)-0.5*ones(6,1)); % update x0
        end
    end
end
