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

d_eps = 0.0000001;
eps_f = 0.00001;
eps_arr   = [0: d_eps : eps_f];
time = linspace(0, 10, size(eps_arr,2));

sig_1 = 0;
eps_p_1 = [0 0 0 0 0 0];
sig_arr(1) = sig_1;
dl = 0;
alpha_1 = 0;
alpha_2 = 0;
ka = 1.0e4;
n = [0 0 0 0 0 0];

for t = 2 : size(eps_arr,2)

  d_eps = eps_arr(t) - eps_arr(t-1);

  eps_2 = [eps_arr(t) 0 0 0 0 0];
  e_2 = eps_2 - (1/3)*(tr(eps_2)) * [1 1 1 0 0 0];
  e_p_1 = eps_p_1 - (1/3)*(tr(eps_p_1)) * [1 1 1 0 0 0];

  eps_p_trial = eps_p_1;
  s_trial = 2*M*(e_2 - e_p_1);
  f_trial = tn(s_trial) - sqrt(2/3) * (Sy + ka * alpha_2);
  dl = 0;

  if (f_trial < 0)
    printf ("LINEAR=YES\n")
    % elastic state
    s = s_trial;
  else
    printf ("LINEAR = NO\n");
    % return mapping
    n = s_trial / tn(s_trial);
    alpha_2 = alpha_1;
    its = 0; 
    do
      g        = - sqrt(2/3) * ka * alpha_2 + tn(s_trial) - 2*M*dl
      dg       = - 2*M;
      dl      -= g/dg;
      alpha_2 += sqrt(2/3) * dl;
      its     += 1;
    until(abs(g)<1.0 || its>4)

  end

  eps_p_2 = eps_p_1 + dl*n;
  sig_2   = K * tr(eps_2) * [1 1 1 0 0 0] + s_trial - 2*M*dl*n;
  eps_p_1 = eps_p_2;
  alpha_1 = alpha_2;

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
