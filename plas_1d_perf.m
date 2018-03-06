% Plasticity for 1D model
% supositions:
% * perfect plasticity

E = 1.0e9;
d_eps = 0.0001;

eps_f = 0.001;
eps_arr   = [0: d_eps : eps_f];
eps_p_arr = zeros(size(eps));
eps_e_arr = zeros(size(eps));
sig_arr   = zeros(size(eps));
time = linspace(0, 10, size(eps_arr,2));

sig_y = 5.0e5;
sig_1 = 0;
eps_p_1 = 0;

for t = 1 : size(eps_arr,2)

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
plot(time, sig_arr,'*-b',"linewidth",2); print -djpg sig.jpg 
