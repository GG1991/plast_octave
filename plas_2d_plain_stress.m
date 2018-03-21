% Test plain stress plasticity 2D model
% > perfect plasticity f = |s|^2 - s_y^2
% > Von Mises criterium

d_eps = 0.0001;
eps_f = 0.001;
%eps_arr = [[0: d_eps : eps_f],[eps_f : -d_eps : -2*eps_f], [-2*eps_f : d_eps : 0]] * [1 0 0; 0 0 0; 0 0 0];
eps_arr = [[0: d_eps : eps_f]; zeros(size([0: d_eps : eps_f])) ; zeros(size([0: d_eps : eps_f]))];

% variables var_x where x is 1 or 2 (old or new) t means "trial" or "test"
eps_e_1 = zeros(3,1); % elastic strain
eps_p_1 = zeros(3,1); % plastic strain

time = linspace(0, 10, size(eps_arr,2));

for t = 2 : size(eps_arr,2)

  [sig_2, eps_e_2, eps_p_2] = func_2d_plain_stress (eps_arr(:,t), eps_e_1, eps_p_1)
  eps_e_1 = eps_e_2;
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
