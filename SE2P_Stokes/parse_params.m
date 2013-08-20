% ------------------------------------------------------------------------------
function p = parse_params(opt)

% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% verify all assumptions on parameters

% 2) grid size the same in x and y
assert(abs(opt.box(2)-opt.box(1))<eps)

% step size
L = opt.box(1);
h = L/opt.M;

% grid point coordinates (x,y)
x = linspace(0,L-h,opt.M);
y = linspace(0,L-h,opt.M);

% Gaussian
P = opt.P;
if( isfield(opt,'m')), m = opt.m; else m = 0.95*sqrt(pi*P); end;
w = h*(P-1)/2;

% domain z, staggered
if isfield(opt,'zLim') 
  % domain in z-direction given 
  a = opt.zLim(1);
  b = opt.zLim(2);
  z = (a+h/2):h:(b-h/2);
  Mz = length(z);
  
  % may not hit z=b exactly, due to the constraint that h_z = h_x = h_y.
  b = z(end)+h/2;
  Lz = b-a;
    
  % is the box big enough?
  assert(a <= -w)
  assert(b >= opt.box(3)+w)
  
else 
  % otherwise: extend box with w to accomodate non-wrapped Gaussians
  z = (-w+h/2):h:(opt.box(3)+w-h/2);
  Lz = opt.box(3) + 2*w;
  Mz = length(z);

  % verify that the z-grid is consistent
  b = z(end)+h/2;
  assert(abs(b-(opt.box(3)+w))<eps) 
end

% sampling factor (oversampling)
if( isfield(opt,'sampling_factor')), s = opt.sampling_factor; else s=4; end;

% algorithm selection
grid_method = 'vanilla';
    
% collect
p.M=opt.M;
p.Mz=Mz;
p.P = P;
p.PH = (P-1)/2;
p.w = w;
p.m = m;
p.box = opt.box;
p.wbox = [0 L; 0 L; -w  -w+Lz];
p.L = L;
p.Lz = Lz;
p.h = h;
p.x = x;
p.y = y;
p.z = z;
p.grid_method = 'vanilla';
p.sampling_factor = s;
