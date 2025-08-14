%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for system ID on particular data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rk, K_err, P_err] = sysid(dim_x, dim_u, x_diffs, xs, us, Q, R, K_true, P_true)
    if rank([xs , us]) < dim_x + dim_u % the rank is deficient, just give up
        rk = 0;
        K_err = 0; 
        P_err = 0;
    else % the rank condition is satisfied, proceed
        rk = 1;
        AB = ([xs , us]\x_diffs).'; % solve the system
        A_id = AB(:,1:dim_x); % separate A
        B_id = AB(:,dim_x+1:end); % separate B
        [K_id,P_id,poles] = lqr(A_id,B_id,Q,R); % solve the LQR problem 
        % assign outputs
        K_err = max(reshape(abs(K_id - K_true),[dim_x*dim_u,1])); % error of K obtained from system ID
        P_err = max(reshape(abs(P_id - P_true),[dim_x*dim_x,1])); % error of P obtained from system ID
    end
end