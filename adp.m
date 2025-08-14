%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method for adp on particular data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rk, K_iterr, P_iterr, K_err, P_err] = adp(dim_x, dim_u, K0, iterations, kron_x_diffs, xxs, xus, Q, R, K_true, P_true)
    if rank([xxs, xus]) < (0.5*dim_x*(dim_x + 1) + dim_x*dim_u) % the rank is deficient, just give up
        rk = 0;
        K_iterr = 0;
        P_iterr = 0; 
        P = zeros(dim_x, dim_x);
        K = zeros(dim_u, dim_x);
        K_err = 0; 
        P_err = 0;
    else % the rank condition is satisfied, proceed
        rk = 1; 
        Ps = zeros(dim_x*dim_x, iterations); % matrix for storing P iterates
        Ks = zeros(dim_x*dim_u, iterations); % matrix for storing K iterates
        Ks(:,1) = reshape(K0,[dim_x*dim_u,1]); % save initial K
        Kk = K0; % initialize K

        for i = 1:iterations 
        
            th = [kron_x_diffs, -2*xxs*kron(eye(dim_x),Kk.'*R) - 2*xus*kron(eye(dim_x), R)]; % build matrix theta
            xi = -xxs*reshape(Q + Kk.'*R*Kk,[dim_x*dim_x,1]); % build matrix xi
            pk = th\xi; % solve system
            P = pk(1:dim_x*dim_x); % separate P
            K = pk(dim_x*dim_x+1:end); % separate K
            Ks(:,i+1) = K; % store the present K
            Ps(:,i) = P; % store the present P
            Kk = reshape(K,[dim_u,dim_x]); % update the present K
        
        end
        % assign outputs
        K_iterr = max(abs(Ks(:,end) - Ks(:,end-1))); 
        P_iterr = max(abs(Ps(:,end) - Ps(:,end-1)));
        P = reshape(Ps(:,end),[dim_x,dim_x]);
        P = 0.5*(P+ P.');
        K = reshape(Ks(:,end),[dim_u,dim_x]);
        K_err = max(reshape(abs(K - K_true),[dim_x*dim_u,1]));
        P_err = max(reshape(abs(P - P_true),[dim_x*dim_x,1]));
    end
end