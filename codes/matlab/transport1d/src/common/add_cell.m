function x = add_cell(a,b)
na = size(a); nta = prod(na); la = length(na);
nb = size(b); ntb = prod(nb); lb = length(nb);
if la ~= lb , error('Array dimensionalities not the same.'); end
nd = nb - na;
if sum(nd) > 0, error('Array sizes not the same.'); end
a = a(:);
b = b(:);
x = cell(nta,1);
for t=1:nta, x{t} = a{t} + b{t}; end
x = reshape(x,na);