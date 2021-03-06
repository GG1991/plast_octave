% Plain stress plasticity for 2D model
% > plasticity with hardening f2 = |s|^2 - (s_y + K*alpha)
% > Von Mises criterium

function [sig_2, eps_e_2, eps_p_2, alpha_2] = func_2d(eps_2, eps_e_1, eps_p_1, alpha_1)

 P       = (1/3)*[2 -1 0 ; -1 2 0 ; 0 0 6];
 E       = 1.0e9;
 nu      = 0.3; 
 G       = E/(2 + 2*nu);
 D0      = (E/(1-nu**2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
 sig_y   = 2.0e11;
 K       = 2.5e20;

 eps_1   = eps_e_1 + eps_p_1;
 d_eps   = eps_2 - eps_1;

 % trial values
 eps_e_t = eps_e_1 + d_eps;
 sig_t   = D0 * eps_e_t;
 eps_p_t = eps_p_1;
 alpha_t = alpha_1;
 S2      = sig_y + K * alpha_1;
 f_2_t   = (1/2) * sig_t' * P * sig_t - (1/3) * S2;

 if (f_2_t <= 0)

   printf("is linear\n");
   eps_e_2 = eps_e_t;
   eps_p_2 = eps_p_t;
   sig_2   = sig_t;
   alpha_2 = alpha_t;

 else

   printf("NOT linear\n");
   A = (1/6)*(sig_t(1) + sig_t(2))^2;
   B = (1/2)*(sig_t(1) - sig_t(2))^2;
   C = 2 * sig_t(3)^2;
   b = 2 * G;

   % begin newton-raphson loop
   dl   = 0.0;
   its  = 0;
   while (its < 50) 
     a    = (1/3) * dl * E / (1-nu);
     S2   = sig_y + K *(alpha_1 + dl); % plasticity with hardening
     phi2 = A/((1+a*dl)^2) + B/((1+b*dl)^2) + C/((1+b*dl)^2);
     q    = (1/2)*phi2 - (1/3)*S2; printf("q = %f\n",q);
     if (abs(q) < 1) break; endif;
     dq   = -(A*a/((1+a*dl)^3) + B*b/((1+b*dl)^3) + C*b/((1+b*dl)^3)) - K/3;
     dl  -= q/dq;
     its += 1;
   endwhile

   sig_2   = inv(inv(D0) + dl*P) * inv(D0) * sig_t;
   eps_p_2 = eps_p_1 + dl * P * sig_2;
   eps_e_2 = eps_2 - eps_p_2;
   alpha_2 = alpha_1 + dl;

 endif

endfunction
