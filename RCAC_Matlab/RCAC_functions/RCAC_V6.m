function [u_out, theta_out] = RCAC_V6(kk, u_in, z_in, yp_in, r_in, FILT, varargin)
%	This function computes vector theta that optimizes the RCAC cost.
%
%   Inputs:     kk      An integer denoting discrete time k.
%               u_in    u(k-1)
%               z_in    z(k-1)
%               yp_in   y(k-1)
%               r_in    r(k-1)
%               FILT    Filter structure containing the following fields
%                       FILT.TYPE   State space or FIR transfer function
%                       FILT.A      A matrix of Gf
%                       FILT.B      B matrix of Gf
%                       FILT.C      C matrix of Gf
%                       FILT.nf     Number of coefficients in the FIR tf
%                       FILT.Nu     FILT.nf lz by lu coefficients stacked
%                                   left to right.
%
%   Outputs:    u_out   u(k) = Phi(k) theta(k)
%               theta   theta(k)
%               OptData Data structure to output anything you want
%
%   Author:     Ankit Goel,
%               Aerospace Engineering Depratment,
%               University of Michigan, Ann Arbor.
%
%   Version:    6.0     2017/11/20
%
%   License:    RCAC license. Need to figure out what this should mean!!


%Persistent varaibles
persistent ltheta lu lz ly         % Vector lengths
persistent u_h yp_h z_h r_h             % Data buffers
% For Batch
persistent P_k theta_k                  % For RLS
persistent intg


% 14 November 2016. ACS2 variables
persistent PHI_window  u_window z_window  pc pn nf Nc
persistent Xphi Xu  % Filtering states
persistent PHI_filt_window  u_filt_window


%%
if (kk == 1)                            % Define lengths
    lu  = size(u_in(:,1),1);
    lz  = size(z_in(:,1),1);
    ly  = size(yp_in(:,1),1);
end

if isempty(varargin)                    % FLAG defaults
    FLAG.RegZ       = 1;                % Regressoin on performance
    FLAG.R0         = 1e+0;             % R_theta = controller gain weight
    FLAG.Rz         = 1;                % Performance weight
    FLAG.Ru         = 0;             % Control weight
    FLAG.lambda     = 1;                % Forgetting factor

    % Controller structure
    FLAG.RegZ       = 1;                % 0->y goes in the controller. 1->z
    FLAG.FF         = 0;                % Feedforward structure
    FLAG.Integrator = 0;                % Integrated performance in the regressor
    FLAG.FIR        = 0;                % FIR structure
    FLAG.ContType   = 'dense';          % Pi is dense

else
    FLAG            = varargin{1};
end

if (kk==1)                              % Controller structure
    Nc      = FLAG.Nc;
    if (isfield(FLAG, 'Nc_ij'))
        Nc      = max(max(max(FLAG.Nc_ij)), Nc);
    end
    nf      = FILT.Nf;
    ltheta  = CalculateRegSize( FLAG.Nc, lu, lz, ly, FLAG);
    u_h     = zeros(lu,Nc);              % Control Buffer
    z_h     = zeros(lz,Nc+1);            % Performance Buffer
    r_h     = zeros(lz,Nc+1);            % Command Buffer
    theta_h = zeros(ltheta,2);          % Controller Parameter Buffer
    intg    = 0;

    if FLAG.RegZ
        yp_h    = zeros(lz,Nc+1);            % Regressor variable Buffer
    else
        yp_h    = zeros(ly,Nc+1);            % Regressor variable Buffer
    end

end






%% Build buffers
u_h(:,2:end)    = u_h(:,1:end-1);       % shift controls by 1
z_h(:,2:end)    = z_h(:,1:end-1);
r_h(:,2:end)    = r_h(:,1:end-1);
yp_h(:,2:end)   = yp_h(:,1:end-1);

u_h(:,1)        = u_in;                 % u(kk-1)
z_h(:,1)        = z_in;                 % z(kk-1)
r_h(:,1)        = r_in;                 % r(kk-1)
if FLAG.RegZ == 1
    yp_h(:,1) = z_in;                   % Performance z in regressor
elseif FLAG.RegZ == 0
    yp_h(:,1) = yp_in;                  % Measurement y in regressor
end


%% Build Regressor

intg    = intg + z_h(:,1);              % g(k) = g(k-1) + z(k-1)
PHI     = BuildRegressor( FLAG.Nc, lu, u_h, yp_h, r_h, intg,FLAG  );
% disp(PHI)
% PHI(k+1) according to aseem

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLLER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (kk == 1)                            % Define lengths
    %% Following variables are used in Gf optimization
    pc          = FLAG.window;
    pn          = round(pc+nf+Nc-1+1);
    PHI_window  = zeros(lu, ltheta, pn);
    u_window    = zeros(lu, pn);
    z_window    = zeros(lz, pc);

    PHI_filt_window     = zeros(lz, ltheta, 2*ltheta);
    u_filt_window       = zeros(lz, 2*ltheta);

    P_k     = eye(ltheta)/FLAG.R0;
    theta_k = 0e-2+zeros(ltheta,1);


end



%% Store data buffers here

nf_end = 5;
u_window(:,2:nf+nf_end)       = u_window(:,1:nf+nf_end-1);
z_window(:,2:nf+nf_end)       = z_window(:,1:nf+nf_end-1);
PHI_window(:,:,2:nf+nf_end)   = PHI_window(:,:,1:nf+nf_end-1);



u_window(:,1)           = u_in;     % the first column vector is u(kk-1)
z_window(:,1)           = z_in;     % the first column vector is z(kk-1)
PHI_window(:,:,1)       = PHI;      % the first lu by ltheta matrix is PHI(kk), used to calculate u(kk)


%%


[ PHI_filt, u_filt, z_filt, Xphi, Xu ] = FilterSignals( kk, FILT, PHI_window, u_window,z_window, Xphi, Xu );

u_filt_window(:,2:end)       = u_filt_window(:,1:end-1);
PHI_filt_window(:,:,2:end)   = PHI_filt_window(:,:,1:end-1);
u_filt_window(:,1)           = u_filt;
PHI_filt_window(:,:,1)       = PHI_filt;






FLAG.PHI = PHI;
if kk> max(ltheta,nf*lz*lu)+1
    control_on 	= 1;

    [ theta_k, P_k ] = RLS_update(kk, theta_k, P_k, PHI_filt, z_filt, u_filt,  FLAG);

else
    control_on  = 0;
end
%     disp([kk z_filt u_filt PHI_filt theta_k'])
%     disp([kk z_filt u_filt])


theta_out   = theta_k;

% Compute the controller
theta_out   = theta_out(:);
u_out       = control_on * PHI * theta_out;    % u(kk) = Phi(kk) theta(kk)



end