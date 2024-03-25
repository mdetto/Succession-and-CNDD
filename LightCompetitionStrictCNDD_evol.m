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



function [n,r,t,R,a,b,lambda] = LightCompetitionStrictCNDD_evol(S,mp,tr)


t0 = 1;
opts = optimset('TolX',1e-18,'TolFun',1e-18);

% lambda = (1-rmin)/(S+1);

n  = zeros(S,1);
t  = zeros(S,1);
r  = zeros(S,1);
R  = zeros(S,1);


%% growth out vs. up
    c = 1.5;
    c1 = c+1;
%     tr = rmin  + rand(S,1)*(1-rmin);
    r(1) = max(tr(tr<1));
%     rp(1) = (r.^-c - 1).^(-1/c);
    t(1) = fminbnd(@fun1,0,t0,opts);
    
%     z = linspace(0.95,t0,1000);
%     x=zeros(1000,1);
%     for i=1:100
%         x(i) = fun1(z(i));
%     end
%     
%     plot(z,x)
    
    R(1) = 1 - t(1);
    n(1) = (t(1)^-c - 1)*r(1)^c1/c;

    i=1;
%     while ks(i)<max(kr) && k(i)<max(kr)
    while t(i)>min(tr) && r(i)>t(i)
        i=i+1;

        r(i) = max(tr(tr<t(i-1) & tr<r(i-1)));
        t(i)  = fminbnd(@fun2,0,t(i-1),opts);
        R(i) = t(i-1) - t(i);
        n(i) = (t(i)^-c - t(i-1)^-c)*r(i)^c1/c;
    end
   
    n=n(1:i);
    r=r(1:i);
    t=t(1:i);
    R=R(1:i);
    
    dt = [1-t(1);-diff(t)];
    a = [1-r(1); t(1:end-1)-r(2:end)];
    b = r-t;


%% plotting
if nargout==0
    
    clf
    subplot(121)
    semilogy(r,n,'o','linewidth',2,'markersize',5);hold all
    plot(r(1),n(1),'ro','markersize',5)
    ylabel('{\itN_i}','FontName','Cambria Math','interpreter','tex')
    xlabel('\it{t}_r','FontName','Cambria Math')
        
    subplot(122)
    plot(r,n./dt,'o','markersize',5)
    yline(1)
    ylim([0.5 1.5])
    
    ylabel('{\itN_i} /  \Delta{\itk_i^*}','FontName','Cambria Math')
    xlabel('\it{t}_r','FontName','Cambria Math')
    
end

%% function for i=1
function y = fun1(t1)

    n1 = (t1^-c - 1)*r(1)^c1/c;
    R1 = n1*(1+mp*n1);
    y = abs(R1 - 1 + t1);

end

%% function for i>1
function y = fun2(ti)
       
    ni = (ti^-c - t(i-1)^-c)*r(i)^c1/c;
    Ri = ni*(1+mp*ni);
    y = abs(Ri - t(i-1) + ti);
end

end

