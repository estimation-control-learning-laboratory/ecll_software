function [ PHI ] = BuildRegressor( Nc, lu, u_h, y_h, r_h, intg_state,FLAG  )
%	This function constructs Phi(k) for various structures.
%
%   Inputs:     Nc      SISO controller order. Not necessarily MIMO
%                       controller order
%               lu      number of controller outputs
%               u_h     [u(k-1) ... u(k-Nc)]
%               y_h     [y(k-1) ... y(k-Nc)]
%               r_h     [r(k-1) ... r(k-Nc)]
%               intg_state  z(k-1)
%               FLAG    specifications for controller structures. Must
%                       contain following fields
%                       FLAG.Nc             An integer
%                       FLAG.Integrator     0 or 1
%                       FLAG.ContType       ManySISO, supersparse, sparse,
%
%   Outputs:    PHI     Phi(k) 
%                       To be used as u(k) = Phi(k) theta(k)
%
%   Author:     Ankit Goel,
%               Aerospace Engineering Depratment, 
%               University of Michigan, Ann Arbor. 
%
%   Version:    6.0     2017/11/20  
%               6.1     2017/12/01      Added ManySISO_Nij for arbitrary 
%                                       controller order in each channel  
%   
%   License:    RCAC license. Need to figure out what this should mean!!


U   = u_h(:,1:Nc);      % U = [u(k-1); ...; u(k-Nc)]
V   = y_h(:,1:Nc);      % V = [xi(k-1); ...; xi(k-Nc)] xi is whatever goes in the regressor
R   = r_h(:,1:Nc);      % R = [r(k-1); ...; r(k-Nc)] used if feedforward is used. Mostly not used

if isfield(FLAG, 'FIR')
    if FLAG.FIR
        U = [];
        V = y_h(:,1:max(1,Nc));
        %FLAG.Integrator = 0;
    end
end


if strcmp(FLAG.ContType, 'dense')
    if FLAG.FF
        phi = [U(:); V(:); R(:)];
    elseif FLAG.Integrator
        phi = [U(:); V(:); intg_state];
    else
        phi = [U(:); V(:)];
    end
    PHI = kron(eye(lu), phi');
elseif strcmp(FLAG.ContType, 'sparse')
    if FLAG.Integrator
        phi_z   = [V(:); intg_state] ;
    else
        phi_z   = [V(:)] ;
    end
    PHI_u   = zeros(lu, lu*Nc);
    for ii = 1:Nc
        phi_u = diag(u_h(:,ii));
        PHI_u(:,(ii-1)*lu+1:(ii-1)*lu+lu) = phi_u;
    end
    PHI_z   = kron(eye(lu), phi_z');
    PHI     = [PHI_u PHI_z];
    if FLAG.FIR
        PHI     = PHI_z;
    end
elseif strcmp(FLAG.ContType, 'supersparse')
    if FLAG.Integrator
        phi_z   = [V(:); intg_state] ;
    else
        phi_z   = [V(:)] ;
    end
    PHI_u   = U;
    PHI_z   = kron(eye(lu), phi_z');
    PHI     = [PHI_u PHI_z];
elseif strcmp(FLAG.ContType, 'PID')
    PHI     = [y_h(:,1) intg_state y_h(:,1)-y_h(:,2)];
elseif strcmp(FLAG.ContType, 'PI')
    PHI     = [y_h(:,1) intg_state];    
elseif strcmp(FLAG.ContType, 'P')
    PHI     = [y_h(:,1)];    
elseif strcmp(FLAG.ContType, 'FSF')
    PHI     = [y_h(:,1)'];
elseif strcmp(FLAG.ContType, 'FSFI')
    PHI     = [y_h(:,1)' -intg_state];
end





end

