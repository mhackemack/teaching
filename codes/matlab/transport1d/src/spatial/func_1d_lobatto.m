function varargout = func_1d_lobatto(varargin)
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
vmin = min(verts(1),verts(2));
vmax = max(verts(1),verts(2));
vdiff = vmax - vmin;
[v,~,~] = lglnodes(order+1);
v = vmin + vdiff*v;

% Compute functions
%-------------------------------------------------------------------------------
if ~gbool
    bout = func_1d_common(v,qx);
else
    [bout,gout] = func_1d_common(v,qx);
end

% Set output arguments
%-------------------------------------------------------------------------------
varargout{1} = bout;
if gbool, varargout{2} = gout; end