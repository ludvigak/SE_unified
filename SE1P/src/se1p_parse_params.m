function popt = se1p_parse_params(opt)
 
% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% copy all params
popt = opt;

% step size
popt.L = opt.box(3);
popt.h = popt.L/popt.M;

% Gaussian
popt.m = .95*sqrt(pi*popt.P);
w = popt.h*popt.P/2;
eta = (2*w*popt.xi/popt.m)^2;
popt.c = 2*popt.xi^2/eta;

% if(eta<1)
%    w = max(w, sqrt((1-eta)/2)*popt.m/popt.xi);
% end

% we want even grids
deltaM = 2*ceil(w/popt.h);
popt.w = popt.h*deltaM/2;
popt.Mx = popt.M+deltaM;
popt.My = popt.M+deltaM;
popt.Lx = popt.box(1)+2*popt.w;
popt.Ly = popt.box(2)+2*popt.w;

% sampling factor (oversampling)
if( isfield(opt,'sg')), popt.sg = opt.sg; else popt.sg=1; end;
if( isfield(opt,'sl')), popt.sl = opt.sl; else popt.sl=1; end;
if( isfield(opt,'s0')), popt.s0 = opt.s0; else popt.s0=1; end;
if( isfield(opt,'nl')), popt.nl = opt.nl; else popt.nl=3; end;
%if( isfield(opt,'R')),  popt.R  = opt.R;  else popt.R=L*sqrt(2); end
if( ~isfield(opt,'k0mod')), popt.k0mod = 1; end% just to define something.

overM = ceil(popt.sg * popt.Mx/2)*2;
popt.sg = overM / popt.Mx;

popt.R = sqrt(popt.Lx^2+popt.Ly^2);

% % find local pad modes to apply higher over sampling factor.
if(popt.sg~=popt.sl)
    n = popt.nl;
    if(n>1)
        n = min(floor((popt.M-1)/2),n); % half modes should be at most
                                        % half of M
    end
    if(mod(popt.M,2)==0)
        popt.local_pad = [2:n+1 popt.M-n+1:popt.M]; % 0 mode is the first element
        popt.k0mod = 1;
    else
        popt.local_pad = [1:n+1 popt.M-n-1:popt.M-1];% 0 mode is the last element
        popt.k0mod = popt.M;
    end
else
    popt.local_pad = 1;
end

% collect
popt.PH = popt.P/2;
wbox = [-popt.w  -popt.w+popt.Lx;-popt.w  -popt.w+popt.Ly; 0 popt.L];
popt.free_offset = wbox(1:2,1);
popt.nl = ceil(numel(popt.local_pad)/2);
popt.eta = (2*popt.w*popt.xi/popt.m)^2;
popt.c = 2*popt.xi^2/popt.eta;