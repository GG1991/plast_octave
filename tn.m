function norm = tn (tensor)

norm = 0;
for i = 1:6
   norm += tensor(i);
endfor

endfunction
