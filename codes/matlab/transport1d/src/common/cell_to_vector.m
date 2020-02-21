function out = cell_to_vector(x)
x = x(:); out = [];
for t=1:length(x), out = [out;x{t}]; end