function out = vector_to_cell(x,stride,dims)
n = dims; nt = prod(n);
out = cell(nt);
for t=1:nt
    s = ((t-1)*stride+1):t*stride;
    out{t} = x(s);
end
out = reshape(out,dims);