% Plasticity for 2D model
% supositions:
% * perfect plasticity f = |s| - s_y 
% * von mises criterium

d_eps = 0.0001;
eps_f = 0.001;
%eps_arr = [[0: d_eps : eps_f],[eps_f : -d_eps : -2*eps_f], [-2*eps_f : d_eps : 0]] * [1 0 0; 0 0 0; 0 0 0];
eps_arr = [0: d_eps : eps_f] * [1 0 0]'; % ex ey 2exy

% variables var_x where x is 1 or 2 (old or new) t means "trial" or "test"
eps_e_1 = zeros(3,1); % elastic strain
eps_e_2 = zeros(3,1);
eps_e_t = zeros(3,1);
eps_p_1 = zeros(3,1); % plastic strain
eps_p_2 = zeros(3,1);
sig_1   = zeros(3,1); % stress tensor
sig_2   = zeros(3,1);
s_1     = zeros(3,1); % deviator s_ij = sig_ij - sig_m * delt_ij
s_2     = zeros(3,1);

P = (1/3)*[2 -1 0 ; -1 2 0 ; 0 0 6];
E = 1.0e9;
nu = 0.3; 
C = (E/(1-nu**2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

time = linspace(0, 10, size(eps_arr,2));

sig_y = 5.0e5;

for t = 2 : size(eps_arr,2)

 d_eps = eps_arr(t) - eps_arr(t-1);
 eps_e_t = eps_e_1 + d_eps;
 sig_t = C * eps_e_t;

end

%figure();
%plot(time, eps_p_arr, '*-r', "linewidth", 2,...
%     time, eps_e_arr, '*-g', "linewidth", 2,...
%     time, eps_arr  , '*-b', "linewidth", 2); print -djpg eps.jpg 
%
%figure();
%plot(eps_arr, sig_arr,'*-b',"linewidth",2); print -djpg sig.jpg 
%
%
%data = [time',eps_arr',eps_e_arr',eps_p_arr',sig_arr'];
%save output.dat -ascii data
