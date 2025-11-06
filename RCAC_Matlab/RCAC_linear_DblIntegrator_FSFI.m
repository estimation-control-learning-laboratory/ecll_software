% clc
clear all
close all
format longe
% format
warning('off','all')
tic
randn('state',3)
% rand('state',2)
addpath './RCAC_functions/'
LoadFigurePrintingProperties

steps           = 2000;


FLAG.Nc         = 4;
FLAG.window     = steps;
FLAG.Rz         = 1;
FLAG.Ru         = 0e+1;

% Controller structure
FLAG.RegZ       = 0;            % 0->y goes in the controller.
FLAG.FF         = 0;
FLAG.Integrator = 1;
FLAG.FIR        = 0;
FLAG.ContType   = 'dense';
FLAG.ContType   = 'PID';
FLAG.ContType   = 'PI';
FLAG.ContType   = 'P';
FLAG.ContType   = 'FSF';
% FLAG.ContType   = 'FSFI';

%% system definition

SysType = 'MP'; 'NMP';
switch SysType
    case 'MP'
        A   = [0.3 0.2; 0.1 0.6];
        B   = [0; 1];
        C   = [1 0];
    case 'NMP'
        [A,B,C,~] = tf2ss(poly(1.2), poly([0.2 -0.5]));
end
Gzu = get_TF_from_ABCD(A,B,C);
% Memory allocation

A = [0 1; 0 0];
B = [0; 1];
C = [1 0];
D = 0;
dt = 0.1
sysC = ss(A,B,C,D)
sysD = c2d(sysC, dt)

A = sysD.A;
B = sysD.B;
C = sysD.C;
D = sysD.D;
[NUM,DEN] = ss2tf(A,B,C,D);
G = tf(NUM,DEN,1);
%%
lx  = size(A,1);
lz  = size(C,1);
ly  = size(C,1);
lu  = size(B,2);
xi  = ones(lu,1);

FLAG.lx = lx;
FLAG.ly = ly;

x0 = randn*ones(lx,1);

R0_array = 10.^[-8:8];
R0_array = 10.^[1:3];
% R0_array = 10.^(8);
% leg_array = {};
for rr = 1:1
    for R0 = R0_array

        %     leg_array{end+1} = num2cell(R0)

        ltheta      = CalculateRegSize( FLAG.Nc, lu, lz, ly, FLAG);
        x           = zeros(lx, steps);
        u           = zeros(lu, steps);
        y0          = zeros(ly,steps);
        z           = zeros(lz,steps);
        theta       = zeros(ltheta,steps);
        r           = zeros(ly,steps);
        intg        = 0;




        % Filter Definition
        %   Gf = N1/q + N2/q^2 + .., where Ni is lz x lu matrix
        %   Gf is specified by FILT.Nu = [N1 N2 ...]

        FILT.TYPE   = 'TF'; 'SS';
        switch SysType
            case 'MP'
                FILT.Nu     = [0 1];
            case 'NMP'
                FILT.Nu     = [1 -1.2];  % Zero polynomial is included
        end
        FILT.Nf     = size(FILT.Nu,2)/lu;

        FLAG.R0     = R0;
        FLAG.lambda = 1;
        FLAG.Ru     = 0e+1;

        r01 = randn;
        r02 = randn;
        r03 = randn;

        %% Simulation
        for ii = 1:steps

            r(:,ii)   = 1 + 1*(ii>steps/2);
            r(:,ii)   = r01 + r02*(ii>steps/3) + r03*(ii>2*steps/3);
            r(:,ii)   = 0;


            if ii == 1
                x(:,ii)     = x0;
                y0(:,ii)    = C*x(:,ii);

                z(:,ii)     = y0(:,ii) - r(:,ii);

                [u(:,ii), theta(:,ii)] = ...
                    RCAC_V6(ii, zeros(lu,1), 0*z(:,ii), x(:,ii), r(:,ii), FILT, FLAG);
            else
                x(:,ii)     = A*x(:,ii-1) + B*u(:,ii-1);
                y0(:,ii)    = C*x(:,ii) + 0e-3*randn;

                z(:,ii)     = y0(:,ii) - r(:,ii);

                [u(:,ii), theta(:,ii)] = ...
                    RCAC_V6(ii, u(:,ii-1), z(:,ii-1), x(:,ii-1), r(:,ii), FILT, FLAG);

            end
            if abs(u(:,ii))>1e3
                break
            end
        end
        [R0 z(:,ii)]
        % theta(:,ii)
        % G_er = 1/(G*theta(:,ii) - 1)
        %%
        time = 1:ii;
        figure(2)
        stairs(time, r(:,1:ii)', 'k--','linewidth',1,'handlevisibility','off')
        hold on; grid on; axis tight
        stairs(time, y0(:,1:ii)', 'linewidth',2)

        %     ylim([-1 3])

        figure(3)
        stairs(time, u(:,1:ii)', 'linewidth',2)
        hold on; grid on; axis tight
        ylabel('$u$')
        %     ylim([-3 3])

        % figure(1)
        %
        % %dTheta = A(up);
        % irow    = 2;
        % icol    = 2;
        % subplot(irow,icol,1)
        %     grid on; axis tight; hold on
        %     stairs(time, (abs(z(:,1:ii)))', 'k', 'linewidth',2)
        %     set(gca,'yscale','log','ytick',10.^[-10:2:10])
        %     ylabel('${\rm{log}}_{10}|z|$')
        %
        %
        % subplot(irow,icol,2)
        %     %stairs(time, log10(abs(u-Au)), 'linewidth',2)
        %     stairs(time, u(:,1:ii)', 'linewidth',2)
        %     grid on; axis tight; hold on
        %     ylabel('$u$')
        %     hold on
        %     ylim([-1 3])
        %
        % subplot(irow,icol,3)
        %     stairs(time, y0(:,1:ii)', 'linewidth',2)
        %     hold on; grid on; axis tight
        %     stairs(time, r(:,1:ii)', 'k','linewidth',2)
        %     ylabel('$y$')
        %     ylim([-1 3])
        % subplot(irow,icol,4)
        %     stairs(time, theta(:,1:ii)', 'linewidth',2)
        %     grid on; axis tight; hold on
        %     ylabel('$\theta$')
        %     ylim([-3 3])
        %ylim([min(theta1(:,ii))-.1 max(theta1(:,ii))+.1  ])



    end
end
%%
figure(2)
% aa = legend( '$10^1$', '$10^2$','$10^3$', '$10^4$','$10^5$', '$10^6$');
% set(aa, 'box','off','location','northwest','interpreter','latex');
% title(['$N_1 = ' num2str(FILT.Nu) ...
%        ', n_{\rm c} = ' num2str(FLAG.Nc)...
%        ', \lambda = ' num2str(FLAG.lambda)...
%        ', R_{u} = ' num2str(FLAG.Ru) '$' ])
title(['$N_1 = ' num2str(FILT.Nu) ...
    ', \lambda = ' num2str(FLAG.lambda)...
    ', R_{u} = ' num2str(FLAG.Ru) '$' ])

ylabel('$y$')
















