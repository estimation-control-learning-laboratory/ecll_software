function [ ltheta ] = CalculateRegSize( Nc, lu, lz, ly,  FLAG)
% This function computes the regressor size for various controller
% structure. Regressor = Phi(k) is lu times ltheta. We compute ltheta in
% this function.



Integrator  = FLAG.Integrator;
ContType    = FLAG.ContType;
if (isfield(FLAG, 'FF'))
    FF = FLAG.FF;
else
    FF = 0;
end
if FF ~= 0
    FF = 1;
end
if Integrator ~= 0
    Integrator = 1;
end

if Nc == 0
    ltheta  = lu * lz;      % Static feedback. Used in Ray's RCMR
else
    if strcmp(ContType, 'dense')
        ltheta  = lu * (lu*Nc + lz*Nc + ly*Nc*FF + lz*Integrator );
    elseif strcmp(ContType, 'sparse')
        ltheta  = ~FLAG.FIR*lu * Nc + lu*lz * Nc + lu*ly*Nc*FF + lu*lz*Integrator;
    elseif strcmp(ContType, 'supersparse')
        ltheta  = Nc + lu*lz * Nc + + lu*ly*Nc*FF + lu*lz*Integrator;
    elseif strcmp(ContType, 'StepInjectionFIR')
        ltheta  = lu*Nc;
    elseif strcmp(ContType, 'ManySISO')
        ltheta  = lu*lz*(2*Nc+Integrator);
    elseif strcmp(ContType, 'ManySISO_Nij')
        if (isfield(FLAG, 'Nc_ij'))
            Nc_ij   = FLAG.Nc_ij;
        else
            Nc_ij   = Nc*ones(lu,lz);
        end
        if (isfield(FLAG, 'Intg_mat'))
            Intg_mat= FLAG.Intg_mat;
        else
            Intg_mat= zeros(lu,lz);
        end
        ltheta  = sum((2*Nc_ij(Nc_ij>0))) + sum(Nc_ij(Nc_ij<0)) + sum(Intg_mat(:));
    elseif strcmp(ContType, 'PID')
        ltheta  = 3;
    elseif strcmp(ContType, 'PI')
        ltheta  = 2;
    elseif strcmp(ContType, 'P')
        ltheta  = 1;
    elseif strcmp(ContType, 'FSF')
        ltheta  = FLAG.lx;
    elseif strcmp(FLAG.ContType, 'FSFI')
        ltheta  = FLAG.lx+FLAG.ly;
    end
end

if FLAG.FIR
    ltheta = lu*lz*Nc;
end


end