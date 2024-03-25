%% LightCompetitionStrict
%
% This code computes the equilibrium density of a light competition model with different tradeoffs. 
% Supplement of 
%  " _Maintenance of high diversity in mechanistic forest dynamics models of
% competition for light_ "
%
% Matteo Detto, Jonathan M. Levine and Stephen W. Pacala
% Department of Ecology and Evolutionary Biology, Princeton University
%


%% Description
% Inputs
%
% * _t0_: fixed disturbace interval
% * _S_ : number of species in the pool
% * options: a structure contain model specific parameters
%
% Outputs
%
% * _k_ : a vector of parameter _k_
% * _n_ : a vector of equilibrium species density
% * _t_ : a vector of closing times
% * _ks_ : a vector of establishment conditions ( _k>ks_ )



function [k,n,ks,t,kk,nk,n0,k0,kmax,omega] = LightCompetitionStrictCNDD_crown(t0,S,mp)


c = 1.5;
% mp = 0.01;
F = 1;
opts = optimset('TolX',1e-18,'TolFun',1e-18);

n  = zeros(S,1);
t  = zeros(S,1);
ks = zeros(S,1);
k  = zeros(S,1);
R  = zeros(S,1);



%% growth out vs. up
    
    c1 = c/(c+1);
    n0 = (sqrt(1+4*mp*F)-1)/(2*mp);
    k0 = 1/(n0*t0^c);
    kmax = k0*10;
%     lambda = kmax./(S+1);
%     k(1) = k0+lambda;
%     t1 = fminbnd(@(x) fun1,0,t0,opts);
%     ks1 = c2/F/t1^c2*t0;
%     S0 = (kmax-ks1)./(2*lambda)+1;
    

    kr = rand(S,1)*kmax;
    k(1) = min(kr(kr>k0));
    t(1) = fminbnd(@fun1,0,t0,opts);
    R(1) = F*c*(log(t0./t(1))) + F;
    n(1) = (sqrt(1+4*mp*R(1))-1)/(2*mp);
    Z = k(1)*n(1);
    ks(1) = max(k(1),1/F/t(1)^c);
    
 
    i=1;n2  = zeros(S,1);
    while ks(i)<max(kr) && ks(i)>=k(i)
        i=i+1;
        k(i) = min(kr(kr>ks(i-1)));
        
        t(i)  = fminbnd(@fun2,0,t(i-1),opts);
        R(i) = F*c*(log(t(i-1)./t(i)));
        n(i) = (sqrt(1+4*mp*R(i))-1)/(2*mp);
        n2(i) = (1-Z.*t(i).^c)./(k(i)*t(i)^c);
        Z = Z+k(i)*n(i);
        ks(i) = max(k(i),1/F/t(i)^c);
        
    end
   
    n=n(1:i);
    k=k(1:i);
    ks=ks(1:i);
    t=t(1:i);
    R=R(1:i);
    
    kk=linspace(k0,kmax,1000);
    nk = F./kk;
    
    dk=diff(ks);
%     plot(k(2:end)-ks(1:end-1),ks(2:end)-k(2:end),'o')
% %     refline(1,0)
%     pause

    

%% plotting
if nargout==0
    
    clf
    subplot(121)
    semilogy(k(1:end)/k0,n(1:end),'o','linewidth',2,'markersize',5)
    hold all;plot(k(1)/k0,n(1),'bo','linewidth',2,'markersize',5)
    ylabel('{\itN_i}','FontName','Cambria Math','interpreter','tex')
    
    subplot(122)

    
    nn = n(2:end)./dk;
    plot(k(2:end)/k0,nn,'o',kk/k0,nk,'r','linewidth',2,'markersize',5)
%     plot(k(2:end)/k0,n(2:end)./dk,'o',kk/k0,nk,'r','linewidth',2,'markersize',5)
%     hold all
%     plot(k(2:end)/k0,R(2:end)./dk,'ro','markersize',5)
    ylabel('{\itN_i} /  \Delta{\itk_i^*}','FontName','Cambria Math')
    xlabel('\it{k / k}_0','FontName','Cambria Math')
    
    

end

%% function for i=2
function y = fun1(t1)
   
    ni = 1./(k(1).*t1.^c);
    Ri = ni + mp*ni.^2;
    y = abs(Ri - F*c*(log(t0./t1)) - F);
    
end

%% function for i>2
function y = fun2(ti)
   
     ni = (1-Z.*ti.^c)./(k(i)*ti^c);
     Ri = ni + mp*ni.^2;
     y = abs(Ri - F*c*(log(t(i-1)./ti)));
end

end
