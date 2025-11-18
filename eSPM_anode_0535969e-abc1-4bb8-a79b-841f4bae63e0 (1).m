clear all
% Applied current density in a galvanostatic setting
global iapp R T FF Npmax ll tau D_particle gamma L1 epsv tplus D_electrolyte
global kappa_el
iapp =2; % A/m^2
%

% Parameters --------------------------------------------------------------
ll = 1.0e-6; % Non-dimensionalization char length
tau = 1; % 
R = 8.314; T =300;
FF = 96500; % Faraday's 
Npmax = 29.4*1.0e3; % Particle occupancy mol/m^3 Measured by experiment
D_particle = 1.0e-11;

kappa_s = 10; %S/m

L1 = 1.0e-4; %m
epsv = 0.5;
tplus = 0.4;
gamma = 1/3;
D_electrolyte = 5.34*1.0e-10; %m^2/s; Measured diffusivity Li+ in electrolyte 
kappa_el = 1.1; %S/m

% -------------------------------------------------------------------------
% Solve for the particle --------------------------------------------------
% -------------------------------------------------------------------------
m = 2;
%
x = linspace(0,10,100); % Radius = 10 microns.
% x is equivalent to r in the notes.
t = linspace(0,2.62*3600/iapp,100000);

% Call the PDE solver of MATLAB and find c(t):
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

c_val = sol(:,end,1); % c at r=R on the radius

%
u_sol = sol(:,:,1);
%
 
V0 = 0.88; % in volts -- to be chosen from expt data
Re = 0.01;
iex = 10;

capacity = iapp*t./3600*1.0e3;

muf = R*T.*log(c_val./(1-c_val));
% -------------------------------------------------------------------------
% Particle solution ends --------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Solve for the electrolyte -----------------------------------------------
% -------------------------------------------------------------------------
m = 0;
%
x = linspace(0,200,200); % length = 200 microns.
% x is equivalent to r in the notes.
t = linspace(0,2.62*3600/iapp,100000);

% Call the PDE solver of MATLAB and find c(t):
sol = pdepe(m,@pdex1pde_el,@pdex1ic_el,@pdex1bc_el,x,t);


%
c_sol = sol(:,:,1);
phi_el = sol(:,1:100,2);
phi_sol = sol(:,:,2);

%
%figure
%plot(x,phi_sol(end,:)*R*T/FF,'b-')
%figure
%plot(x,c_sol(end,:),'r-')



% Electrolyte solution ends -----------------------------------------------
% -------------------------------------------------------------------------


% Calculate cell voltage --------------------------------------------------
% Numerically integrate to find phi_av
phi_av = [];
for i=1:1:length(t)
    phi_av(i)=(ll/L1)*trapz(x(1:100),phi_el(i,:)); 
end

% eSPM Vcell
V_cell = V0-muf./FF-iapp*(L1/3/kappa_s)+phi_av'*R*T/FF+2*R*T/FF*asinh(-iapp/iex/2); 

% SPM Vcell
V_cell_SPM = V0-muf./FF+2*R*T/FF*asinh(-iapp/iex/2); 



figure
plot(capacity,V_cell,'g-')
xlabel('Capacity (mA-hr)','Interpreter','latex');
ylabel('Half cell voltage (V)','Interpreter','latex');
hold on
plot(capacity,V_cell_SPM,'g--')










% Particle subroutines ----------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
global iapp R T FF Npmax ll tau
%

iapp_nd = iapp*tau/(ll*Npmax*FF); % iapp_hat in the notes
% 
pl = 0.0; % Is not important when m=2
ql = 0.0; % Is not important when m=2
pr = -iapp_nd;
qr = 1.0;
end

function u0 = pdex1ic(x)
u0 = 0.001; % initial condition on c
end

function [c,f,s] = pdex1pde(x,t,u,DuDx)
%
global iapp R T FF Npmax ll tau D_particle
%
% Diffusivity of Li in graphite particle
%
Df = (D_particle/ll^2)*tau; % Non-dimensionalized diffusivity
%
%
Diffus = Df; % Parameter
c = 1;
f = Diffus/(1-u)*DuDx;
s = 0;
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc_el(xl,ul,xr,ur,t)
%
global iapp R T FF Npmax ll tau
%
%
iapp_nd = iapp*tau/(ll*Npmax*FF); % iapp_hat in the notes
% 
pl = [0.0; 0.0]; 
ql = [1.0; 1.0]; 
pr = [-iapp_nd; ur(2)];
qr = [1.0; 0.0];
end


% -------------------------------------------------------------------------

function u0 = pdex1ic_el(x)
u0 = [0.5; 0.0]; % initial condition on c
end

% -------------------------------------------------------------------------

function [c,f,s] = pdex1pde_el(x,t,u,DuDx)
global iapp R T FF Npmax ll tau gamma L1 epsv tplus D_electrolyte kappa_el
%
%
iapp_nd = iapp*tau/(ll*Npmax*FF); % iapp_hat in the notes
%
%
kappa_hat = kappa_el*tau*R*T/(Npmax*FF^2*ll^2);
Df = (D_electrolyte/ll^2)*tau; % Non-dimensionalized diffusivity

Diffus = Df; % Parameter

if x<=100
c = [epsv*gamma ; 0];    
f = [Diffus*DuDx(1)+tplus*iapp_nd*x*(ll/L1); ...
     kappa_hat*(DuDx(2)-2*(1-tplus)*(1/u(1))*DuDx(1))];
s = [-(ll/L1)*iapp_nd; -iapp_nd];
else
c = [gamma ; 0];    
f = [Diffus*DuDx(1)+tplus*iapp_nd; ...
     kappa_hat*(DuDx(2)-2*(1-tplus)*(1/u(1))*DuDx(1))];
s = [0; 0];
end
end


