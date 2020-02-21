function varargout = func_1d_common(pts, x)
varargout{1} = get_1D_values(pts, x);
if nargout > 1, varargout{2} = get_1D_gradients(pts, x); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_1D_values(v, x)
ord = length(v)-1;
if ord == 0
    out = ones(length(x),1);
elseif ord == 1
    out = [(v(2)-x), (x-v(1))]/(v(2) - v(1));
elseif ord == 2
    out = [(x-v(2)).*(x-v(3))/(v(1)-v(2))/(v(1)-v(3)),...
           (x-v(1)).*(x-v(3))/(v(2)-v(1))/(v(2)-v(3)),...
           (x-v(1)).*(x-v(2))/(v(3)-v(1))/(v(3)-v(2))];
elseif ord == 3
    out = [(x-v(2)).*(x-v(3)).*(x-v(4))/(v(1)-v(2))/(v(1)-v(3))/(v(1)-v(4)),...
           (x-v(1)).*(x-v(3)).*(x-v(4))/(v(2)-v(1))/(v(2)-v(3))/(v(2)-v(4)),...
           (x-v(1)).*(x-v(2)).*(x-v(4))/(v(3)-v(1))/(v(3)-v(2))/(v(3)-v(4)),...
           (x-v(1)).*(x-v(2)).*(x-v(3))/(v(4)-v(1))/(v(4)-v(2))/(v(4)-v(3))];
else
    out = ones(length(x),ord+1);
    for i=1:ord+1
        for j=1:ord+1
            if i==j, continue; end
            out(:,i) = out(:,i).*(x-v(j))/(v(i)-v(j));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_1D_gradients(v, x)
ord = length(v)-1;
if ord == 0
    out = zeros(length(x),1);
elseif ord == 1
    out = ones(length(x),1) * [-1, 1]/(v(2) - v(1));
elseif ord == 2
    out = [  (x - v(2))/((v(1) - v(2))*(v(1) - v(3))) + (x - v(3))/((v(1) - v(2))*(v(1) - v(3))), ...
           - (x - v(1))/((v(1) - v(2))*(v(2) - v(3))) - (x - v(3))/((v(1) - v(2))*(v(2) - v(3))), ...
             (x - v(1))/((v(1) - v(3))*(v(2) - v(3))) + (x - v(2))/((v(1) - v(3))*(v(2) - v(3)))];
elseif ord == 3
    v1 = v(1); v2 = v(2); v3 = v(3); v4 = v(4);
    out = [   ((v2 - x).*(v3 - x))/((v1 - v2)*(v1 - v3)*(v1 - v4)) + ((v2 - x).*(v4 - x))/((v1 - v2)*(v1 - v3)*(v1 - v4)) + ((v3 - x).*(v4 - x))/((v1 - v2)*(v1 - v3)*(v1 - v4)),...
            - ((v1 - x).*(v3 - x))/((v1 - v2)*(v2 - v3)*(v2 - v4)) - ((v1 - x).*(v4 - x))/((v1 - v2)*(v2 - v3)*(v2 - v4)) - ((v3 - x).*(v4 - x))/((v1 - v2)*(v2 - v3)*(v2 - v4)),...
              ((v1 - x).*(v2 - x))/((v1 - v3)*(v2 - v3)*(v3 - v4)) + ((v1 - x).*(v4 - x))/((v1 - v3)*(v2 - v3)*(v3 - v4)) + ((v2 - x).*(v4 - x))/((v1 - v3)*(v2 - v3)*(v3 - v4)),...
            - ((v1 - x).*(v2 - x))/((v1 - v4)*(v2 - v4)*(v3 - v4)) - ((v1 - x).*(v3 - x))/((v1 - v4)*(v2 - v4)*(v3 - v4)) - ((v2 - x).*(v3 - x))/((v1 - v4)*(v2 - v4)*(v3 - v4))];
else
    nx = length(x); klist = (1:ord+1)';
    out = zeros(nx,ord+1); tx = ones(nx,1);
    for i=1:ord+1
        tlist = klist; tlist(i) = []; denom = 1;
        for k=1:ord
            ttx = tx;
            for kk=1:ord
                if k==kk, continue; end
                ttx = ttx.*(x-v(tlist(kk)));
            end
            out(:,i) = out(:,i) + ttx;
        end
        % Calculate denominator
        for j=1:ord
            jj = tlist(j);
            denom = denom*(v(i) - v(jj));
        end
        out(:,i) = out(:,i) / denom;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%