% Plasticity for 2D model
% supositions:
% * perfect plasticity f = |s|^2 - s_y^2
% * von mises criterium

d_eps = 0.0001;
eps_f = 0.001;
%eps_arr = [[0: d_eps : eps_f],[eps_f : -d_eps : -2*eps_f], [-2*eps_f : d_eps : 0]] * [1 0 0; 0 0 0; 0 0 0];
eps_arr = [[0: d_eps : eps_f]; zeros(size([0: d_eps : eps_f])) ; zeros(size([0: d_eps : eps_f]))];

% variables var_x where x is 1 or 2 (old or new) t means "trial" or "test"
eps_e_1 = zeros(3,1); % elastic strain
eps_e_2 = zeros(3,1);
eps_e_t = zeros(3,1);
eps_p_1 = zeros(3,1); % plastic strain
eps_p_2 = zeros(3,1);
sig_1   = zeros(3,1); % stress tensor
sig_2   = zeros(3,1);
sig_t   = zeros(3,1);
s_1     = zeros(3,1); % deviator s_ij = sig_ij - sig_m * delt_ij
s_2     = zeros(3,1);

P = (1/3)*[2 -1 0 ; -1 2 0 ; 0 0 6];
E = 1.0e9;
nu = 0.3; 
G = E/(2 + 2*nu);
D0 = (E/(1-nu**2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

time = linspace(0, 10, size(eps_arr,2));

sig_y = 5.0e11;

for t = 2 : size(eps_arr,2)

 d_eps   = eps_arr(:,t) - eps_arr(:,t-1);
 eps_e_t = eps_e_1 + d_eps;
 sig_t   = D0 * eps_e_t;
 eps_p_t = eps_p_1 ;
 S2      = sig_y; % perfect plasticity
 f_2_t   = (1/2) * sig_t' * P * sig_t - (1/3) * S2;

 if (f_2_t <= 0)
   printf("is linear\n");
   eps_e_2 = eps_e_t;
   sig_2   = sig_t;
   f_2     = f_2_t;
 else
   printf("is NOT linear\n");
   A = (1/6)*(sig_t(1) + sig_t(2))**2;
   B = (1/2)*(sig_t(1) - sig_t(2))**2;
   C = 2 * sig_t(3)**2;
   b    = 2*G;

   % begin newton-raphson loop
   dl   = 0.0;
   its  = 0;
   while (its < 20) 
     a    = (1/3) * dl * E / (1-nu);
     S2   = sig_y; % perfect plasticity
     phi2 = A/((1+a*dl)**2) + B/((1+b*dl)**2) + C/((1+b*dl)**2);
     q    = (1/2)*phi2 - (1/3)*S2;
     printf("q = %f\n",q);
     if (q < 0.001) break; endif;

     dq   = -(1/2)*(2*A*a/((1+a*dl)**3) + 2*B*b/((1+b*dl)**3) + 2*C*b/((1+b*dl)**3));
     dl  -= q/dq;

     its += 1;
   endwhile

 endif

 eps_e_1 = eps_e_2;
 sig_1 = sig_2;

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
