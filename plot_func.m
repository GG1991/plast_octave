
P       = (1/3)*[2 -1 0 ; -1 2 0 ; 0 0 6];
E       = 1.0e9;
nu      = 0.3; 
G       = E/(2 + 2*nu);
D0      = (E/(1-nu**2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
sig_y   = 2.0e11;
factor  = 1.0e8

eps_e_1 = [0 0 0]';
eps_p_1 = [0 0 0]';
eps_2   = [0.005 0 0]';
eps_1   = eps_e_1 + eps_p_1;
d_eps   = eps_2 - eps_1;

% trial values lets see what happen
eps_e_t = eps_e_1 + d_eps;
sig_t   = D0 * eps_e_t;
eps_p_t = eps_p_1;

S2      = sig_y; % perfect plasticity

A = (1/6)*(sig_t(1) + sig_t(2))**2;
B = (1/2)*(sig_t(1) - sig_t(2))**2;
C = 2 * sig_t(3)**2;
b = 2 * G;

% begin newton-raphson loop
its  = 0;
dl0 = 0;
dlf = 0.001;
dl = linspace(dl0,dlf,1000);
q  = zeros(size(dl));

for i = 1:size(dl,2);
  a    = (1/3) * dl(i) * E / (1-nu);
  S2   = sig_y; % perfect plasticity
  phi2 = A/((1+a*dl(i))**2) + B/((1+b*dl(i))**2) + C/((1+b*dl(i))**2);
  q(i) = ((1/2)*phi2 - (1/3)*S2)/factor;
  if (abs(q) < 0.00001) break; endif;
  dq   = -((A*a/((1+a*dl(i))^3) + B*b/((1+b*dl(i))^3) + C*b/((1+b*dl(i))^3)))/factor;
endfor

%-----

max_its = 8;
q_tri = zeros(1,max_its);
dl_tri = zeros(1,max_its);
dl_tri(1) = 0.0001;
a    = (1/3) * dl_tri(1) * E / (1-nu);
S2   = sig_y; % perfect plasticity
phi2 = A/((1+a*dl_tri(1))^2) + B/((1+b*dl_tri(1))^2) + C/((1+b*dl_tri(1))^2);
q_tri(1) = ((1/2)*phi2 - (1/3)*S2)/factor;
dq   = -(A*a/((1+a*dl_tri(1))^3) + B*b/((1+b*dl_tri(1))^3) + C*b/((1+b*dl_tri(1))^3))/factor;

for i = 2 : max_its
  dl_tri(i) = dl_tri(i-1) - q_tri(i-1)/dq;
  a    = (1/3) * dl_tri(i) * E / (1-nu);
  S2   = sig_y; % perfect plasticity
  phi2 = A/((1+a*dl_tri(i))^2) + B/((1+b*dl_tri(i))^2) + C/((1+b*dl(i))^2);
  q_tri(i) = ((1/2)*phi2 - (1/3)*S2)/factor;
  dq   = -(A*a/((1+a*dl_tri(i))^3) + B*b/((1+b*dl_tri(i))^3) + C*b/((1+b*dl_tri(i))^3))/factor;
endfor

fig=figure();
xl = [dl0 dlf];
yl = [0   0  ];  
grid on;
plot(dl_tri, q_tri, '-g', "linewidth", 1); print -djpg q.jpg 
%plot(xl, yl, '-r', dl, q, '-b', "linewidth", 2); print -djpg q.jpg 
data = [dl_tri', q_tri'];
save q_dl_tri.dat -ascii data
data = [dl', q'];
save q_dl.dat -ascii data
