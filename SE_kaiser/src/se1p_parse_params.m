function popt = se1p_parse_params(opt)
 
% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% copy all params
popt = opt;

% step size
popt.L = opt.box(1);
popt.h = popt.L/popt.M;

% Gaussian
popt.m = .95*sqrt(pi*popt.P);
popt.w = popt.h*popt.P/2;
popt.eta = (2*popt.w*popt.xi/popt.m)^2;
popt.c = 2*popt.xi^2/popt.eta;
popt.p_half = (mod(opt.P,2)==0)*opt.P/2+(mod(opt.P,2)~=0)*(opt.P-1)/2;

popt.My = ceil(popt.box(2)/popt.h)+popt.P;
popt.Mz = ceil(popt.box(3)/popt.h)+popt.P;
popt.Ly = popt.box(2)+2*popt.w;
popt.Lz = popt.box(3)+2*popt.w;

% Even grids if even grids are in the z direction
if(mod(popt.My,2)~=0)
    popt.My = 2*ceil(popt.My/2);
end
if(mod(popt.Mz,2)~=0)
    popt.Mz = 2*ceil(popt.Mz/2);
end
    
popt.Ly = popt.h*popt.My;
popt.Lz = popt.h*popt.Mz;

% h should be the same in all directions
assert(abs(popt.Ly/popt.My-popt.h)<eps)
assert(abs(popt.Lz/popt.Mz-popt.h)<eps)

% sampling factor (oversampling)
if( isfield(opt,'s')), popt.s = opt.s; else popt.s=1; end;
if( isfield(opt,'s0')), popt.s0 = opt.s0; else popt.s0=2; end;
if( isfield(opt,'n')),
    popt.n = min(opt.n,ceil(popt.M/2));
else 
    popt.n=max(ceil(popt.M/2),1); 
end;

if(isfield(opt,'beta'))
    popt.beta=opt.beta*popt.P; 
else 
    popt.beta=2.3*popt.P; 
end

popt.R = sqrt(popt.Ly^2+popt.Lz^2);

% increase s and such that FFTN has integer size vectors.
popt.s  = max(ceil(popt.s *popt.My)/popt.My,ceil(popt.s *popt.Mz)/popt.Mz);
popt.s0 = max(ceil(popt.s0*popt.My)/popt.My,ceil(popt.s0*popt.Mz)/popt.Mz);


% local pads
% FIXME: We assume that the same number of modes in x and y directions
% are oversampled.
n = popt.n;
if(n>1)
    n = min(floor((popt.M-1)/2),n); % half modes should be at most
                                    % half of M
end

popt.local_pad = [2:n+1 popt.M-n+1:popt.M];

% collect
popt.PH = popt.P/2;
wbox = [0 popt.L; -popt.w -popt.w+popt.Ly; -popt.w  -popt.w+popt.Lz];
popt.free_offset = wbox(2:3,1);
