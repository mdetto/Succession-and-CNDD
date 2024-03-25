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
% * _JM_: number of species in the metacommunity
% * _alfa_ : CNDD parameter
% * _rmin_ : lower lint for break-even time

%
% Outputs
%
% * _k_ : a vector of parameter _k_
% * _n_ : a vector of equilibrium species density
% * _t_ : a vector of closing times
% * _r_ : a vector of break-even times



function [n,r,t,R,a,b,lambda] = LightCompetitionStrictCNDD_gost(JM,alfa,rmin)



opts = optimset('TolX',1e-18,'TolFun',1e-18);

% lambda = (1-rmin)/(JM+1);

% n  = zeros(JM,1);
t  = zeros(JM,1);
r  = zeros(JM,1);
R  = zeros(JM,1);

lambda = JM/(1-rmin);
%% growth out vs. up
    c = 1.5;
    c1 = c+1;
%     tr = rmin  + rand(JM,1)*(1-rmin);
%     r(1) = max(tr(tr<1));
r(1) = 1+log(rand)/lambda;
%     rp(1) = (r.^-c - 1).^(-1/c);
t(1) = fminbnd(@fun1,0,r(1),opts);
    
%     z = linspace(0.95,1,1000);
%     x=zeros(1000,1);
%     for i=1:100
%         x(i) = fun1(z(i));
%     end
%     
%     plot(z,x)
    
    R(1) = 1 - t(1);
%     n(1) = R(1)/(1+alfa*R(1));

    i=1;
%     while ks(i)<max(kr) && k(i)<max(kr)
%     while t(i)>min(tr) && r(i)>t(i)
    while t(i)>rmin
        i=i+1;

%         r(i) = max(tr(tr<t(i-1) & tr<r(i-1)));
        r(i) = t(i-1)+log(rand)/lambda;
        if r(i)>rmin
        t(i)  = fminbnd(@fun2,0,r(i),opts);
%         if r(i)<t(i)
%             z = linspace(0,t(i-1),10000);
%             x = zeros(size(z));
%             for j=1:length(z)
%                 x(j) = fun3(z(j));
%             end
% 
%             plot(log(z),x)
% 
%         end
        R(i) = t(i-1) - t(i);
%         n(i) = R(i)/(1+alfa*R(i));
        else
            t(i)=0;
        end
    end
    use = r>rmin;
    r=r(use);
    t=t(use);
    R=R(use);
    
    n = R./(1+alfa*R);
    a =  [1-r(1); t(1:end-1)-r(2:end)];
    b = r-t;


%% plotting
if nargout==0
    
    dt = R;
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
    R1 = n1/(1-alfa*n1);
    y = abs(R1 - 1 + t1);

end

%% function for i>1
function y = fun2(ti)
       
    ni = (ti^-c - t(i-1)^-c)*r(i)^c1/c;
    Ri = ni/(1-alfa*ni);
    y = abs(Ri - t(i-1) + ti);
end


end

