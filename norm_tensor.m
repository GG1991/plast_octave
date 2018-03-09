function [norm] = norm_tensor(tensor)

norm = 0.0;
for i = 1 : size(tensor,1)
  for j = 1 : size(tensor,2)
    norm += tensor(i,j)*tensor(i,j);
  end
end

endfunction
