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



function [r,n,t,index] = LightCompetitionStrictCNDD_linear(tr,mp)

tr = tr(tr>=0);
S = length(tr);
n  = zeros(S,1);
t  = zeros(S,1);
r  = zeros(S,1);
c1 = 2.5;

r(1) = max(tr(tr<1));
% n(1) = 2*(1-r(1));
n(1) = 2*c1/(c1+2*mp*r(1))*(1-r(1));
t(1) = 1 - n(1);

i=1;
while t(i)>min(tr) && r(i)>t(i)
    i=i+1;

    r(i) = max(tr(tr<t(i-1) & tr<r(i-1)));
%     n(i) = 2*(t(i-1)-r(i));
    n(i) = 2*c1/(c1+2*mp*r(i))*(t(i-1)-r(i));
    t(i) =  t(i-1) - n(i);

end

n=n(1:i);
r=r(1:i);
t=t(1:i);
index = ismember(tr,r);

    



