FLAG.RegInit    = 0; % Compute regressor before all data is available
FLAG.FiltInit   = 1; % Compute filter before all data is available
FLAG.ContInit   = 0; % Compute controller before all data is available
FLAG.k_0        = 2+nstr;                  %ANKIT: What is this??

% Controller Type
FLAG.uFIR       = 0;
FLAG.uZERO      = 0;
FLAG.Prop       = 1;
FLAG.Int        = 0;
FLAG.IntProp    = 0;

% Cost function type
FLAG.JInst      = 0;

% Regression Variable
FLAG.RegZ       = 1;

% Forgetting Factor


% Compute regressor vector size
if isfield(FLAG, 'RegZ') && FLAG.RegZ
    lv = lz;
else
    lv = ly;
end


%% Compute controller order and size of control parameter
if FLAG.uZERO
    ltheta  = lu*ly;
    lphi    = ly;
elseif FLAG.uFIR
    if FLAG.Prop
        ltheta  = lu*(Nc+1)*ly;
        lphi    = ly*(Nc+1);
    else
        ltheta  = Nc*(lu*ly);
        lphi    = ly*Nc;
    end
else        % For Proper (instead of strictly proper) controller
    if FLAG.Prop
        ltheta  = lu*(Nc*lu +(Nc+1)*ly);
        lphi    = lu*Nc+ly*(Nc+1);
    else
        ltheta  = Nc*(lu*(ly+lu));
        lphi    = (lu+ly)*Nc;
    end
end

if FLAG.Int % Integrator Size
    ltheta = ltheta+lu*lz;
    lphi = lphi+lz;
end
ltheta  = ltheta+1;
lphi    = lphi+1;

% % FLAG.Rz         = 1*eye(lz);            % Performance penalty
% FLAG.Rf         = 0*eye(lz);            % Filtered Control penalty
% FLAG.Ru         = 0*eye(lu);            % Control penalty
% FLAG.Rtheta     = 1e+2*eye(ltheta);     % Transient Penalty
FLAG.Rdelta     = 0*eye(ltheta);     % Delta theta Penalty
FLAG.theta_0    = 0.0*ones(ltheta,1);
% FLAG.lambda     = 0.999;                    % Must be <= 1







