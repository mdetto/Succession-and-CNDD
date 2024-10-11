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

% dN/dt = R(1 - alfa R) - m N

function [n,r,t,R,a,b,lambda] = LightCompetitionStrictCNDD_R(S,mp,rmin)


t0 = 1;
opts = optimset('TolX',1e-18,'TolFun',1e-18);

lambda = (1-rmin)/(S+1);

n  = zeros(S,1);
t  = zeros(S,1);
r  = zeros(S,1);
R  = zeros(S,1);


%% growth out vs. up
    c = 1.5;
    c1 = c+1;
    tr = rmin  + rand(S,1)*(1-rmin);
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
        n(i) = (t(i)^-c - t(i-1)^-c)*r(i)^c1/c;
    end
   
    n=n(1:i);
    r=r(1:i);
    t=t(1:i);
    R=R(1:i);
    
    dt = [1-t(1);-diff(t)];
    a =  [1-r(1); t(1:end-1)-r(2:end)];
    b = r-t;


%% function for i=1
function y = fun1(t1)


    n1 = (1-t1)*(1-mp*(1-t1));
    y = abs((t1^-c - 1)*r(1)^c1/c - n1);
    
end

%% function for i>1
function y = fun2(ti)
       
     ni = (t(i-1)-ti)*(1-mp*(t(i-1)-ti));
     y = abs((ti^-c - t(i-1)^-c)*r(i)^c1/c - ni);
end



end

