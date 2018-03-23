% Test plain stress plasticity 2D model
% > perfect plasticity f = |s|^2 - s_y^2
% > Von Mises criterium

d_eps = 0.00001;
eps_f = 0.0008;
eps_arr = [[0: d_eps : eps_f]; 0*[0: d_eps : eps_f] ; 0*[0: d_eps : eps_f]];

% variables var_x where x is 1 or 2 (old or new) t means "trial" or "test"
eps_e_1 = zeros(3,1); % elastic strain
eps_e_2 = zeros(3,1);
eps_p_1 = zeros(3,1); % plastic strain
eps_p_2 = zeros(3,1);
sig_2   = zeros(3,1);

time = linspace(0, 10, size(eps_arr,2));
sig_arr_esc = zeros(size(time));
eps_arr_esc = zeros(size(time));
eps_arr_esc(1) = 0; 
sig_arr_esc(1) = 0;

for t = 2 : size(eps_arr,2)

  [sig_2, eps_e_2, eps_p_2] = func_2d_plain_stress (eps_arr(:,t), eps_e_1, eps_p_1);
  eps_e_1 = eps_e_2;
  eps_p_1 = eps_p_2;
  eps_e_arr_esc(t) = norm(eps_e_2);
  eps_p_arr_esc(t) = norm(eps_p_2);
  eps_arr_esc(t)   = norm(eps_arr(:,t));
  sig_arr_esc(t)   = norm(sig_2);

end

figure();
plot(time, eps_p_arr_esc, '*-r', "linewidth", 2,...
     time, eps_e_arr_esc, '*-g', "linewidth", 2,...
     time, eps_arr_esc  , '*-b', "linewidth", 2); print -djpg eps.jpg 
%
figure();
plot(eps_arr_esc, sig_arr_esc,'*-b',"linewidth",2); print -djpg sig.jpg 
%
%
%data = [time',eps_arr',eps_e_arr',eps_p_arr',sig_arr'];
%save output.dat -ascii data
