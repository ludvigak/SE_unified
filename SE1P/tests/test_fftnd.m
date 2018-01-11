clear

sl = 3;
s0 = sl;
local_pad = 1:8;
k0mod = 1;
rep = 8;
Mlist = [16 32 48 64 80 96 112 128];

[t1,t2,e1,e2] = deal([]);

for M = Mlist; 
    H = rand(M,M,M);
    tic
    for k=1:rep
        [G, Gres, G0]= fftnd(H, M, M, M, sl, s0, 1, local_pad, k0mod, 3);
        F = ifftnd(G,Gres,G0, M, M, M, sl, s0, 1, local_pad, k0mod,3);
    end
    t1(end+1)=toc/rep;
    e1(end+1)=max(max(max(abs(F(:)-H(:)))));
    
    tic
    for k=1:rep
        G = fftn(H, [M*sl M*sl M]);
        F = ifftn( G );
        F = F(1:M,1:M,1:M);
    end
    t2(end+1) = toc/rep;
    e2(end+1)=max(max(max(abs(F(:)-H(:)))));
    table(t1',t2',t2'./t1',e1',e2','VariableNames',{'fftnd','fftn','speedup','e1','e2'})
end
