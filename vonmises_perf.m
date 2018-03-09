% Von-mises perfect plasticity

function [sig_2, eps_p_2] = vonmises_perf(eps_1, eps_2, sig_1, eps_p_1)

G = 1.0e10;
sig_y = 1.0e6;
eps_dev_1 = eps_1 - trace(eps_1)/3*eye(3);
eps_dev_2 = eps_2 - trace(eps_2)/3*eye(3);
sig_dev_1 = sig_1 - trace(sig_1)/3*eye(3);

sig_dev_tr = sig_dev_1 + 2*G*(eps_dev_2 - eps_dev_1);

f_tr = norm_tensor(sig_dev_tr) - sqrt(2/3)*sig_y;

if (f_tr < 0)
  sig_dev_2 = sig_dev_tr;
  printf("material is linear\n");
else
  printf("material is in plastic zone\n");
end

sig_2 = sig_dev_2 + trace(sig_dev_2)*3*eye(3);

endfunction
