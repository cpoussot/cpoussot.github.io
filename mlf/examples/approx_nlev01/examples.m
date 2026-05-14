function [Hf,Phi,ip,p_c,p_r,C,sol,XLIM,YLIM] = examples(NUM)
    
    XLIM    = [-1 1];
    YLIM    = [-1 1];
    sol     = @(x) [];
    switch NUM
        case 1
            % True known solution
            sol = @(x) [x; sqrt(1-x); -sqrt(1-x)];
            % Coutour
            C   = .6*exp(1i*linspace(0,2*pi,1e3));
            %
            E   = eye(3);
            A   = @(x) [0 1 0; 1-x(:,1) 0 0; 0 1 x(:,1)];
            Phi = @(x) x(:,1)*E-A(x(:,2));
            H   = @(x) ones(1,3)*(Phi(x)\ones(3,1));
            Hf  = @(x1,x2) H([x1,x2]);
            % IP
            ip{1} = .6*exp(1i*linspace(0,pi,22)); % then complex conjugated
            ip{1} = ip{1}(2:end-1);
            %ip{1} = .6*exp(1i*linspace(1e-2,pi-1e-2,20)); % then complex conjugated
            ip{2} = linspace(.75,1.25,40);
            %ip{2} = linspace(1.25,1.5,40);
            % interlace
            for ii = 1:numel(ip)
                p_c{ii} = ip{ii}(2:2:end);
                p_r{ii} = ip{ii}(1:2:end);
            end
            % complex conjugate ip{1}
            pc{1} = [];
            pr{1} = [];
            for ii = 1:length(p_c{1})
                pc{1} = [pc{1} p_c{1}(ii) conj(p_c{1}(ii))];
                pr{1} = [pr{1} p_r{1}(ii) conj(p_r{1}(ii))];
            end
            p_c{1}  = pc{1};
            p_r{1}  = pr{1};
            ip{1}   = [p_c{1} p_r{1}];
            %
            XLIM    = [-.8 1.6];
            YLIM    = .8*[-1 1];
        case 2
            % Coutour
            C   = .075*exp(1i*linspace(0,2*pi,1e3));
            % 
            N   = 10;
            I   = eye(N);
            E   = diag(logspace(-4,10,N));
            d   = @(x) x(:,1)+0.01*exp(-x(:,1).*x(:,2));
            Phi = @(x) ( d(x).*I+E );
            H   = @(x) ones(1,N)*(Phi(x)\ones(N,1) );
            Hf  = @(x1,x2) H([x1,x2]);
            % IP
            ip{1} = .1*exp(1i*linspace(0,pi,22)); % then complex conjugated
            ip{1} = ip{1}(2:end-1);
            %ip{1} = .1*exp(1i*linspace(.01,pi-.01,20)); % then complex conjugated
            ip{2} = linspace(30,35,40);
            %ip{2} = linspace(20,50,40);
            ip{2} = 18:52;%
            % interlace
            for ii = 1:numel(ip)
                p_c{ii} = ip{ii}(2:2:end);
                p_r{ii} = ip{ii}(1:2:end);
            end
            % complex conjugate ip{1}
            pc{1} = [];
            pr{1} = [];
            for ii = 1:length(p_c{1})
                pc{1} = [pc{1} p_c{1}(ii) conj(p_c{1}(ii))];
                pr{1} = [pr{1} p_r{1}(ii) conj(p_r{1}(ii))];
            end
            p_c{1}  = [pc{1}];
            p_r{1}  = [pr{1}];
            ip{1}   = [p_c{1} p_r{1}];
            %
            XLIM    = .15*[-1 1];
            YLIM    = .15*[-1 1];
        case 3 
            % Coutour
            %Eli = @(z,c,sr,si) (real(z-c)./sr).^2 + (imag(z-c)./si).^2;
            %C   = Eli(1i*linspace(0,2*pi,1e3),-3,2.5,10);
            %[X,Y,Z] = ellipsoid(-3,0,0,2.5,10,0,1e3);
            C   = Ell(-3,2.5,10,linspace(0,2*pi,1e3));
            % 
            zh  = @(x) sqrt( x(:,1).^2 + 2.*x(:,1).*x(:,2) );
            Phi = @(x) [-sinh(x(:,1)./4)         sinh(zh(x)./4)            cosh(zh(x)./4)            0; ...
                        -x(:,1).*cosh(x(:,1)./4) zh(x).*cosh(zh(x))        zh(x).*sinh(zh(x))        0; ...
                        0                        -sinh(3.*zh(x)./4)        -cosh(3.*zh(x)./4)        sinh(x(:,1)./4); ...
                        0                        -zh(x).*cosh(3.*zh(x)./4) -zh(x).*sinh(3.*zh(x)./4) -x(:,1).*cosh(x(:,1)./4)];
            H   = @(x) (ones(1,4))*( Phi(x)\ones(4,1) );
            Hf  = @(x1,x2) H([x1,x2]);
            % IP
            ip{1} = Ell(-3,3,11,linspace(0,pi,22)); % then complex conjugated
            ip{1} = ip{1}(2:end-1);
            ip{2} = linspace(3,4,25);
            % interlace
            for ii = 1:numel(ip)
                p_c{ii} = ip{ii}(2:2:end);
                p_r{ii} = ip{ii}(1:2:end);
            end
            % complex conjugate ip{1}
            pc{1} = [];
            pr{1} = [];
            for ii = 1:length(p_c{1})
                pc{1} = [pc{1} p_c{1}(ii) conj(p_c{1}(ii))];
                pr{1} = [pr{1} p_r{1}(ii) conj(p_r{1}(ii))];
            end
            p_c{1}  = pc{1};
            p_r{1}  = pr{1};
            ip{1}   = [p_c{1} p_r{1}];
            %
            %XLIM    = [-6 0];
            XLIM    = [-6 6];
            YLIM    = [-15 15];
    end
end

%%%
function C = Ell(c,a,b,t)
    %t   = linspace(0,2*pi,N);
    x   = a*cos(t);
    y   = b*sin(t);
    C   = (x+c)+1i*y;
end