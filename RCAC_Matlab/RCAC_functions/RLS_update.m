function [ theta_k, P_k] = RLS_update( kk, theta_k, P_k, PHI_filt, z_filt, u_filt,  FLAG)

%   This function computes the theta and P updates for several forgetting
%   schemes.
%   Default is standard forgetting, where all singular values of the
%   covariance matrix P are divided by the forgetting factor lambda.
%   FLAG.RLS_SVDF   = 0 is Standard Forgetting
%                   = 1 divides ltheta-rank(S) smallest singular values.
%                   = 2 divides ltheta-rank(S) singular values with the
%                       largest drop
%                   = 3 is exponential exponential forgetting factor
%                   = 4 is linear exponential forgetting factor
%                   = 5 divides all but saturates the singular value
%                   = 6 Partitions singular values and minimizes a metric
%                   = 7 is exponential resetting and forgetting from
%                       Salgado's paper
%                   = 8 is my brainchild. Muhahahahaha!!!



lz      = size(z_filt,1);
PHI     = FLAG.PHI;
lu      = 1;
ltheta  = size(theta_k,1);

if ~isfield(FLAG, 'Q')
    FLAG.Q = 0;             % Kalman filter modification to RLS
end

if FLAG.Ru==0
    %     P_k         = P_k - (PHI_filt * P_k )' * inv(FLAG.lambda*eye(lz)+PHI_filt * P_k * PHI_filt') ...
    %         *  (PHI_filt * P_k) + FLAG.Q*eye(ltheta);
    P_k         = P_k - (PHI_filt * P_k )' *  (PHI_filt * P_k) / (1+PHI_filt * P_k * PHI_filt');
%     norm(P_k-P_k')
%     P_k = 0.5*(P_k+P_k');
    
    %P_k     = P_k/FLAG.lambda;
    
    
    
    theta_k = theta_k - P_k * PHI_filt' * (PHI_filt*theta_k + z_filt - u_filt);
    
    
else
    
    Rbarinv = blkdiag(1/FLAG.Rz*eye(lz), 1/FLAG.Ru*eye(lu));
    Rbar    = blkdiag(FLAG.Rz*eye(lz), FLAG.Ru*eye(lu));
    P_k     = P_k - P_k * [PHI_filt; PHI]' * ...
        pinv(FLAG.lambda*Rbarinv +  [PHI_filt; PHI] * P_k * [PHI_filt; PHI]') *...
        [PHI_filt; PHI] * P_k;
    P_k     = P_k/FLAG.lambda;
    theta_k = theta_k - P_k * [PHI_filt; PHI]' * Rbar * [(PHI_filt*theta_k + z_filt - u_filt); PHI*theta_k];
    
end

end

