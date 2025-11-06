function [ PHI_filt, u_filt, z_filt, Xphi, Xu ] = FilterSignals( kk, FILT, PHI_window, u_window, z_window, Xphi, Xu )
%	This function filters the regressor and past inputs used by RCAC
%
%   Inputs:     FILT    Filter structure containing the following fields
%                       FILT.TYPE   State space or FIR transfer function
%                       FILT.A      A matrix of Gf
%                       FILT.B      B matrix of Gf
%                       FILT.C      C matrix of Gf
%                       FILT.nf     Number of coefficients in the FIR tf
%                       FILT.Nu     FILT.nf lz by lu coefficients stacked
%                                   left to right.
%               PHI_window  lu by ltheta by (nf + RD(Gf) + 1) regressor matrix
%               u_window    lu by (nf + RD(Gf) + 1) past controls matrix
%               z_window    lz by (nf + RD(Gf) + 1) past performance matrix
%               Xphi    State of the system driven by Phi(k)
%               Xu      State of the system driven by u(k)
%
%   Outputs:    PHI_filt    Phi_f(k)
%               u_filt      u_f(k)
%               z_filt      z_f(k)
%               Xphi    State of the system driven by Phi(k)
%               Xu      State of the system driven by u(k)
%
%   Author:     Ankit Goel,
%               Aerospace Engineering Depratment,
%               University of Michigan, Ann Arbor.
%
%   Version:    1.0     2017/11/20
%
%   Comments:
%               FILT.A is n  x n
%               FILT.B is n  x lu
%               FILT.C is lz x n
%               the first column vector of u_window is u(kk-1)
%               the first column vector of z_window is z(kk-1)
%               the first lu by ltheta matrix of PHI_window is PHI(kk), used to calculate u(kk)
%
persistent PHI_filt_window u_filt_window
nf      = FILT.Nf;
lz      = size(z_window,1);
ltheta  = size(PHI_window,2);
if kk==1
    PHI_filt_window     = zeros(lz, ltheta, nf);
    u_filt_window       = zeros(lz, nf);
end


GfRD        = 0;            % 1 picks up z(kk-1), Ub(kk-1) and Phi_b(kk-1)
Phi_b_rr    = collapse(PHI_window(:,:,1+GfRD+1:1+nf+GfRD));
U_b_rr      = vec( u_window(:,1+GfRD:1+nf+GfRD-1 ) );
Phi_f_b_rr  = collapse(PHI_filt_window(:,:,1:nf));
U_b_f_rr    = vec( u_filt_window(:,1:nf) );

if ~isfield(FILT, 'Du')
    FILT.Du  = zeros(lz,lz*nf);
end

PHI_filt    = FILT.Nu * Phi_b_rr - FILT.Du * Phi_f_b_rr;
u_filt      = FILT.Nu * U_b_rr - FILT.Du * U_b_f_rr;

PHI_filt_window(:,:,2:end) = PHI_filt_window(:,:,1:end-1);
u_filt_window(:,2:end) = u_filt_window(:,1:end-1);

PHI_filt_window(:,:,1) = PHI_filt;
u_filt_window(:,1) = u_filt;

z_filt      = z_window(:,1+0*GfRD);

% disp([kk z_filt u_filt PHI_filt])


end

