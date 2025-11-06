clc
clear all
close all
format long
warning('off','all')
tic
addpath './RCAC_functions/'


steps           = 10000;
FLAG.steps      = steps;

FLAG.Nc         = 2;
FLAG.window     = steps; 300;

FLAG.optimType  = 'RLS';
FLAG.Rz         = 1;
FLAG.Ru         = 0e+1;

% Controller structure
FLAG.RegZ       = 1;            % 0->y goes in the controller.
FLAG.FF         = 0;
FLAG.Integrator = 1;
FLAG.FIR        = 0;
FLAG.ContType   = 'dense';
%% Van Der Pol system parameters


lx = 2;
lu = 1;
lz = 1;
ly = 1;

xi          = ones(lu,1);
ltheta      = CalculateRegSize( FLAG.Nc, lu, lz, ly, FLAG);
x           = zeros(lx, steps);
u           = zeros(lu, steps);
y0          = zeros(ly,steps);
z           = zeros(lz,steps);
theta       = zeros(ltheta,steps);
r           = zeros(ly,steps);
intg        = 0;

%% Controller settings (FLAG)

FILT.TYPE   = 'TF'; 'SS';
FILT.Nu     = 1;
FILT.Nf     = size(FILT.Nu,2)/lu;

FLAG.R0     = 1e-2;
FLAG.lambda = 1;


%% Simulation

r   = zeros(ly,steps)+0.8;

x0  = [1 1]';
Cd  = [1 1];
for ii = 1:steps
%     r(ii) = 2.5e-4;

    if ii == 1
        %% Initialization
        x(:,1) = VanDerPol(x0, 0,1);
        u(:,1) = zeros(lu,1);
        y(:,1) = Cd*x(:,1);
        z(:,1) = y(:,1) - r(:,ii);
        
        [u(:,ii), theta(:,ii)] = ...
            RCAC_V6(ii, zeros(lu,1), 0*z(:,ii), 0*z(:,ii), r(:,ii), FILT, FLAG);
        
    else
        x(:,ii) = VanDerPol(x(:,ii-1), u(:,ii-1),1);
        y(:,ii) = Cd*x(:,ii);
        z(:,ii) = y(:,ii) - r(:,ii);
        
        [u(:,ii), theta(:,ii)] = ...
            RCAC_V6(ii, u(:,ii-1), z(:,ii-1), 0*z(:,ii-1), r(:,ii), FILT, FLAG);
        
    end

    u(:,ii) = Saturate(u(:,ii), -1, 2);
end

%%
close all
tic

time = 1:steps;


subplot(2,2,1)
    semilogy(time, abs(z))
    grid on; axis tight
    ylabel('z(k)')
    xlabel('Time (seconds)')
    
subplot(2,2,4)
    plot(time, u)
    grid on; axis tight
    ylabel('u(k)')
    xlabel('Time (seconds)')
    
subplot(2,2,3)
    plot(time, y, time, r, 'k', 'linewidth',2)
    grid on; axis tight
    ylabel('y(k)')
    xlabel('Time (seconds)')
    
    
subplot(2,2,2)
    stairs(time, theta')
    grid on; axis tight
    ylabel('\theta(k)')
    xlabel('Time (seconds)')
    
    
    





figure
stairs(x(1,:), x(2,:), 'b') 
grid on
hold on
plot(x(1,1), x(2,1), 'rd', 'linewidth', 4,  'markersize', 10) 
plot(x(1,end), x(2,end), 'bd', 'linewidth', 4,  'markersize', 10) 
xlabel('x_1')
ylabel('x_2')

