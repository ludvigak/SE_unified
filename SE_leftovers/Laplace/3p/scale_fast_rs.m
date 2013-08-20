clear all, close all

xi = 8;

test = '2';

switch test
    case '1'
        N = 1000;
        L = 1;
        S = 5;
        
        for i = 1:20
            [x f] = SE_state(N(i),[L(i) L(i) L(i)]);
            %tic
            [u t(i)] = fast_rs_sum(x, f, xi, L(i), S(i), true);
            %t(i) = toc;
            S(i+1) = S(i)+1;
            L(i+1) = L(i) + L(1)/S(1);
            N(i+1) = round(L(i+1)^3*N(1)/L(1)^3);
        end
        S = S(1:end-1);
        L = L(1:end-1);
        N = N(1:end-1);
        
        plot(N,t,'-*')
        
    case '2'
        
        NN = [500 1500];
        Lbk = [3 2];
        sty = {'-+','-*'}
        
        figure(); hold on;
        for k=1:length(NN)
        N = NN(k);
        L = 1;
        S = 4; % rc = 1/S
        Lb = Lbk(k)*L;
        clear t;
        for i = 1:(1+S(1)*(Lb/L-1))
            [x f] = SE_state(N(i),[L(i) L(i) L(i)]);
            %tic
            [u t(i)] = fast_rs_sum(x, f, xi, L(i), S(i), true);
            %t(i) = toc;
            S(i+1) = S(i)+1;
            L(i+1) = L(i) + L(1)/S(1);
            N(i+1) = round(L(i+1)^3*N(1)/L(1)^3);
        end
        S = S(1:end-1);
        L = L(1:end-1);
        N = N(1:end-1);
        fprintf('%d \t%.2f\n',[N' L']')
        %plot(N,[t.sum]+([t.bx]+[t.verlet])/10,sty{k})
        plot(N,[t.sum],sty{k})
        end
        
        publication_fig
        xlabel('N')
        ylabel('time (s)')
        grid on;
        %fname = 'output/scale_fast_rs';
        %write_fig(1,fname);
end