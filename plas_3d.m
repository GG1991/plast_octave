% Plasticity for 3D model
% supositions:
% * plasticity with hardening f = |s| - s_y 
% * von mises criterium

E  = 1.0e9;
nu = 0.3;
L  = (nu*E)/((1+nu)*(1-2*nu));
K  = E/(3*(1-2*nu));
M  = (3/2) * (K - L);
Sy = 5.0e5;

d_eps = 0.0001;
eps_f = 0.001;
eps_arr   = [0: d_eps : eps_f];
time = linspace(0, 10, size(eps_arr,2));

sig_1 = 0;
eps_p_1 = [0 0 0 0 0 0];
sig_arr(1) = sig_1;

for t = 2 : size(eps_arr,2)

  d_eps = eps_arr(t) - eps_arr(t-1);

  eps_2 = [eps_arr(t) 0 0 0 0 0];
  e_2 = eps_2 - (1/3)*(tr(eps_2)) * [1 1 1 0 0 0];

  eps_p_trial = eps_p_1;
  s_trial = 2*M*(eps_2 - eps_p_1);
  f_trial = tn(s_trial) - sqrt(2/3) * (Sy);

  if (f_trial < 0)
    % elastic state
    printf ("LINEAR=YES\n")
    s = s_trial;
    eps_p_2 = eps_p_trial;
  else
    % return mapping
    printf ("LINEAR = NO\n");
  end

  %sig_2   = sig_2;
  eps_p_1 = eps_p_2;

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
