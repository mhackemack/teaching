function varargout = func_1d_bernstein(varargin)
% Grab input arguments
%-------------------------------------------------------------------------------
verts = varargin{1};
order = varargin{2};
qx    = varargin{3};
gbool = (nargout > 1);

% Perform some error checking
%-------------------------------------------------------------------------------
nv = length(verts);
if nv ~= 2, error('1D cell requires 2 vertices.'); end
if order < 1, error('FEM order must be >= 1.'); end

% Grab geometric information
%-------------------------------------------------------------------------------
vmin  = min(verts(1),verts(2));
vmax  = max(verts(1),verts(2));
vdiff = vmax - vmin;
qref  = (qx-vmin)/vdiff;
t     = qref;
t1    = 1-qref;

% Compute basis functions
%-------------------------------------------------------------------------------
bout = zeros(length(qx),order+1);
for i=0:order
    bout(:,i+1) = binomial(i,order)*(t.^i).*(t1.^(order-i));
end

% Compute gradients
%-------------------------------------------------------------------------------
if gbool
    n    = order;
    gout = zeros(length(qx),n+1);
    gout(:,1)   = -n*t1.^(n-1);
    gout(:,end) = n*t.^(n-1);
    for i=1:n-1
        gout(:,i+1) = binomial(i,n)*(i*t.^(i-1).*(t1).^(n-i) - (n-i)*(t.^i).*(t1.^(n-i-1)));
    end
end

% Set output arguments
%-------------------------------------------------------------------------------
varargout{1} = bout;
if gbool, varargout{2} = gout; end