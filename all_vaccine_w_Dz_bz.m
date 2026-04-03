% Red Queen vaccine effect on evolution speed and infection incidence. All vaccinations, width of peak w is self-consistent.
% Last vaccine delay Dz, vaccine interval Dz1, fixed bz, no restrictions on w, Dx, bz

%% Parameters
R0=1.8;             % basic reproduction number without vaccine, H3N2
ez=0.4;             % vaccine max. efficacy
py=0.024;           % susceptivility to homologous strain
by=4.1;  %1.7;  % K(u) scale, best-fit for different g(s)
%Dz1=1.3;  % Vaccination interval in antigenic coordinate
%% choose the fixed parameter(s) and adjust plot labels in the end
bz=1.7;
Dz=1;
%%
Ub=1e-3; %(-4.5);       % ben. mutation rate, best fit for different g(s)
N=1e8;
U=Ub;               % total mutation rate in asymmetry \beta
distribution_s='uniform'; %'constant';%'halfgaussian';

w20=2;            % initial value of w^2  
du=0.1;    
u=0:du:3000; % antigenic coordinates (denoted -u) for natural cross-immunity K 
u2=-3000:du:3000; % a.c. reflected for virus cross-immunity Kz


% natural immunity crossimmunity (norm. transmission rate)
K=(1+(1-py)/py*exp(-u/by)).^(-1);   % eq. 8
Kder=(1-py)/py/by*exp(-u/by).*(1+(1-py)/py*exp(-u/by)).^(-2);  
cz=0:0.01:1; %  coverage 

%% loop in vaccine breadth b_z

%% choose the parameter that varies between curves
%Yfor Dz=[1 0.5 0 -1]
%for bz=[0.8 1.7 4.1]
for Dz1=[0.32 0.65 1.3]
    
%% Vaccine crossimmunity and its derivative: added shift by zz, made symmetric/asymmetric
    
    zz=(-Dz-Dz1*(0:40))';        % coordinates of 40 vaccines starting from the most recent with respect to virus peak
    U = ones(length(zz),1)*u2 - zz*ones(1,length(u2));      % matrix of shifted coordinates
    Kz = (1+ez/(1-ez)*exp(-abs(U)/bz)).^(-1);  % eq. 10
    Kzder = ez/(1-ez)/bz *(2*(U > 0)-1).*exp(-abs(U)/bz).*(1+ez/(1-ez)*exp(-abs(U)/bz)).^(-2);  

%% Iterating in w
    w2=w20*ones(1,length(cz)); % initial
    w2old=10*w2;
    Re=ones(1,length(cz)); logRe_der=Re; % taking space
    iw=0;
    while max(abs(w2old./w2-1)) > 1e-4 && iw < 100
        iw=iw+1;
% Reproduction number: average over u, product over z
        for k=1:length(cz)
            Omega=1-cz(k)*sqrt(1/(2*pi*w2(k)))*du*exp(-u2.^2/2/w2(k))*(1-Kz)';
            Re(k)=R0*prod(Omega);   
        end
%% Normalization of the recovered,
        AA=zeros(size(cz));
        sigma1=zeros(size(cz));
        for k=1:length(cz)
% Recovered individual density  
            Anew=1; A=77; i3=0;
            while abs(A/Anew-1)>1e-4 && i3<200
                i3=i3+1; A=Anew;
                r=exp(-A*Re(k)*du*cumsum(K));  % eq. S5, denoted s
                totalr=sum(r)*du; 
                Anew=1/totalr;
            end
            r=A*r;
            AA(k)=A; 
% Selection coefficient from natural infections but with reduced Re
            sigma1(k)=Re(k)*du*sum(Kder.*r); % eq. S7, first term
        end % loop in cz
        sigma0=sigma1(1);   % sigma without vaccine, cz=0
        A0=AA(1);           % normalization factor without vaccine
        
%% Finding full sigma including vaccination: average over u, product over z
        for k=1:length(cz)
           logRe_der(k)=sum(cz(k)*sqrt(1/(2*pi*w2(k)))*du*exp(-u2.^2/2/w2(k))*Kzder'./Omega);  
        end
        sigma=sigma1 + logRe_der; 

%% Update w2
        w2old=w2;
% Initial interation of speed and fraction of infected finf
        V=2*sigma*log(N); 
        finf=V.*AA;
% Iterating finf
        fold=0.1*finf; 
        i_f=0;
        while max(abs(finf/fold-1)) > 1e-3 && i_f<100
            i_f=i_f+1;
            fold=finf;
            Ninf=N*finf;
            switch distribution_s
                case 'constant'
                    for i1=1:20 % iterating in V
                        V=2*sigma.*log(Ninf./sqrt(V.^2 .*log(V/Ub)/(sigma.^3*Ub)))./(log(V/exp(1)/Ub).^2+1); % if V > sigma
                    end
                    w2=V./sigma; % variance in antigenic coordinate
                    s_ast=sigma;
                case 'exponential' 
                    v=2*sigma.^2 .*log(Ninf*Ub);     % Initial iteration of adaptation speed
                    s_ast = sqrt(2*v.*log(2/Ub*sqrt(v/2/pi)));  % Most probable fixed allele selection coefficient
                    for i1=1:10
                        xc = s_ast+v./sigma;
                        s_ast = sqrt(2*v.*log(2/Ub*sqrt(v/2/pi)/(1+sigma./xc+v./sigma./s_ast)));
                        v=2*sigma.*(-s_ast+sigma.*log(N*Ub*xc.^2 ./v-1+2*xc.*sigma./v+2*sigma.^2 ./v));
                    end
                    V=v./s_ast;
                    w2=v./s_ast.^2;    % variance in antigenic coordinate
                case 'halfgaussian'
                    s_ast=sigma*sqrt(pi).*log(sigma/Ub).^0.5;   % Initial iteration
                    kappa=log(2*Ninf.*s_ast*Ub/pi./sigma);
                    for i1=1:10
                        kappa=log(2*Ninf.*s_ast*Ub.*(1+kappa).^1.5/pi./sigma./kappa)./...
                        log(s_ast/Ub.*sqrt(kappa./(1+kappa))); 
                % rederived on 17.12.2024 for half-Gaussian rho(s) with average s=1 above
                        s_ast=sigma.*sqrt(pi*kappa./(1+kappa)).*log(s_ast/Ub.*sqrt(kappa./(1+kappa))).^0.5;  % rederived on 17.12.2024
                    end
                    v=pi*sigma.^2/2 .*kappa;
                    V=v./s_ast;
                    w2=v./s_ast.^2;    % variance in antigenic coordinate
                case 'uniform'
                    alpha=10; % step steepness
                    s0=sigma*2*gamma(1+1/alpha)/gamma(1+2/alpha); % to have mean(s)=sigma
                    X=sqrt(alpha)*s0/Ub;
                    % adaptation rate
                    v=2*s0.^2 .*log(Ninf.*s0.*sqrt(log(Ninf.*s0)))./log(X.*log(X).^0.5).^2 ./(1+log(X)/4 ./log(Ninf.*s0)).^2;
                    % effective selection coefficient
                    s_ast=s0.*(sqrt(2)/alpha*log(s0/Ub)).^(1/(alpha-1));
                    V=v./s_ast;
                    w2=v./s_ast.^2;    % variance in antigenic coordinate
            end  % switch distribution s
            finf=V.*AA; % updating finf
        end % iterations in finf
    end % iterations in w2
    
%% Plotting results for positive effective coefficient
ii=find(sigma>0 & [1, sigma(1:(end-1))]>0 & Re > 1);
%ii=find(finf>0 & [1, finf(2:end-1)]>0);
figure(1)
subplot(2,2,1)
    plot(cz(ii), V(ii)/V(1))
    %% choose parameter on curves bz or Dz
    text(cz(ii(end)),V(ii(end))/V(1),sprintf('%g', Dz1))
    hold on
    xlabel('Coverage c_z')
    ylabel('Substitution rate')
    %% choose fixed parameter bz or Dz
    title(sprintf('R0=%g, e_z=%g, p_y=%g, b_y=%g, bz=%g, Dz=%g\n N=%g, U_b=%g, %s, Dz1 on curves' ,R0,ez,py,by,bz,Dz,N,Ub,distribution_s ))

subplot(2,2,2)
    plot(cz(ii), finf(ii)/finf(1))
     %% choose parameter on curves bz or Dz
    text(cz(ii(end)),finf(ii(end))/finf(1),sprintf('%g', Dz1))
    hold on
    xlabel('Coverage c_z')
    ylabel('Incidence')
subplot(2,2,3)
    plot(cz(ii), Re(ii))
    text(cz(ii(end)),Re(ii(end)),sprintf('%g', Dz1))
        hold on
    xlabel('Coverage c_z')
    ylabel('Reproduction number R_e')
end % loop in Dz
figure(1)
subplot(2,2,1)
hold off
box off
subplot(2,2,2)
hold off
box off
subplot(2,2,3)
hold off
box off




    

