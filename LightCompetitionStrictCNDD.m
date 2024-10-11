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



function [k,n,ks,t,R,n0,k0,ks1] = LightCompetitionStrictCNDD(t0,S,mp,kmax)


F = 1;
opts = optimset('TolX',1e-18,'TolFun',1e-18);

n  = zeros(S,1);
t  = zeros(S,1);
ks = zeros(S,1);
k  = zeros(S,1);
R  = zeros(S,1);


%% growth out vs. up
    c = 1.5;
    c1 = c/(c+1);
    c2 = c+1;
    
    if mp ==0 
        n0 = F/c2;
    else
    n0 = (sqrt(1+4*mp*F/c2)-1)/(2*mp);
    end
    k0 = 1/(n0*t0^c);
%     kmax = k0*2;

    lambda = kmax./(S+1);
    k(1) = k0+lambda;
    t1 = fminbnd(@fun1,0,t0,opts);
    ks1 = c2/F/t1^c2*t0;
    
    

%     kr = rand(S,1)*kmax;
    tmin = (c2./kmax).^(1/c2);
    tr = tmin  + rand(S,1)*(1-tmin);
    kr = c2./tr.^c2;
    k(1) = min(kr(kr>k0));
    t(1) = fminbnd(@fun1,0,t0,opts);
    R(1) = F*c1*(t0 - t(1))/t0 + F/c2;
%     n(1) = (sqrt(1+4*mp*R(1))-1)/(2*mp);
    n(1) = 1./(k(1).*t(1).^c);
    Z = k(1)*n(1);
    ks(1) = c2/F/t(1)^c2*t0;
    
    i=1;
    while ks(i)<max(kr) && k(i)<max(kr)
        i=i+1;
        k(i) = min(kr(kr>ks(i-1) & kr>k(i-1)));
        
        t(i)  = fminbnd(@fun2,0,t(i-1),opts);
        R(i) = F*c1*(t(i-1) - t(i))/t0;
%         n(i) = (sqrt(1+4*mp*R(i))-1)/(2*mp);
        n(i) = (1-Z.*t(i).^c)./(k(i)*t(i)^c);
        Z = Z+k(i)*n(i);
        ks(i) = c2/F/t(i)^c2*t0;
    end
   
    n=n(1:i);
    k=k(1:i);
    ks=ks(1:i);
    t=t(1:i);
    R=R(1:i);
    
    kk=linspace(k0,kmax,1000);
%     nk = c1*(n0/t0)^c1*kk.^(-(c+2)/(c+1));
    nk = c*c2^(-(2*c+1)/c2)*(F/t0)^c1*kk.^(-(c+2)/c2);
    
    
    dk=diff(ks);
    
%     a = k(2:end)-ks(1:end-1);
%     b = ks(2:end)-k(2:end);
%     
%     y = t(1:end-1);
%     z = k(2:end);
% %     n2 = (y.^(c+1)./(c+1)-1./z)./(mp./z.^2+1/2*y.^(2*c+1));
%     
%     x1 = 3*c/2*y.^-c.*(1-sqrt(8/3*(2*c+1)./z.*y.^(-c-1) - (13*c+5)/(3*c+3)))/(2*c+1);
%     x2 = 2*c./(c+1)^2.*y.*(z-ks(1:end-1));
%     plot(n(2:end),x2./k(2:end),'o');refline(1,0)
%  
% %     plot(k(2:end)-ks(1:end-1),ks(2:end)-k(2:end),'o')
% %     refline(1,0)
% % 
% % z1=n(2:end)./dk;
% % z2 = c1*(n0/t0)^c1*k(2:end).^(-(c+2)/(c+1));
% % H = mean(z1./z2);
% 
% H = -mean(log(c1*(n0/t0)^c1*k(2:end).^(-(c+2)/(c+1))./n(2:end).*dk)/log(1+mp)); 

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

    n1 = 1./(k(1).*t1.^c);
    R1 = n1 + mp*n1.^2;
    y = abs(R1 - F*c1*(t0-t1)/t0 - F/c2);

end

%% function for i>2
function y = fun2(ti)
   
     ni = (1-Z.*ti.^c)./(k(i)*ti^c);
     Ri = ni + mp*ni.^2;
     y = abs(Ri - F*c1*(t(i-1)-ti)/t0);
end

end

