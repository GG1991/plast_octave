% Plasticity for 3D model
% supositions:
% * perfect plasticity f = |s| - s_y 
% * von mises criterium

E = 1.0e9;
d_eps = 0.0001;

eps_f = 0.001;
%eps_arr   = [[0: d_eps : eps_f],[eps_f : -d_eps : -2*eps_f], [-2*eps_f : d_eps : 0]] * [1 0 0; 0 0 0; 0 0 0];
eps_arr   = [[0: d_eps : eps_f]] * [1 0 0; 0 0 0; 0 0 0];
eps_p_arr = zeros(size(eps));
eps_e_arr = zeros(size(eps));
sig_arr   = zeros(size(eps));
time = linspace(0, 10, size(eps_arr,2));

sig_y = 5.0e5;
sig_1 = 0;
eps_p_1 = 0;
sig_arr(1) = sig_1;
eps_p_arr(1) = eps_p_1;
eps_e_arr(1) = eps_arr(1) - eps_p_1;

for t = 2 : size(eps_arr,2)

  d_eps = eps_arr(t) - eps_arr(t-1);
  sig_trial   = sig_1 + E*d_eps;
  eps_p_trial = eps_p_1;
  f_trial = abs(sig_trial) - sig_y;

  if (f_trial < 0)
    % elastic state
    sig_2   = sig_trial;
    eps_p_2 = eps_p_trial;
  else
    % return mapping
    d_gamma = f_trial/E;
    sig_2   = sig_trial   - d_gamma*E*sign(sig_trial);
    eps_p_2 = eps_p_trial + d_gamma*sign(sig_trial);
  end
  sig_1   = sig_2;
  eps_p_1 = eps_p_2;

  sig_arr(t) = sig_2;
  eps_p_arr(t) = eps_p_2;
  eps_e_arr(t) = eps_arr(t) - eps_p_2;
  
end

figure();
plot(time, eps_p_arr, '*-r', "linewidth", 2,...
     time, eps_e_arr, '*-g', "linewidth", 2,...
     time, eps_arr  , '*-b', "linewidth", 2); print -djpg eps.jpg 

figure();
plot(eps_arr, sig_arr,'*-b',"linewidth",2); print -djpg sig.jpg 


data = [time',eps_arr',eps_e_arr',eps_p_arr',sig_arr'];
save output.dat -ascii data
