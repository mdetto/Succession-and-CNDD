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



function [r,n,t] = LightCompetitionStrictCNDD_linear_v2(S,mp,rmin)

% n  = zeros(S,1);
t  = zeros(ceil(S),1);
r  = zeros(ceil(S),1);

c1 = 2.5;
a = 1/2;
b = mp/c1;

csi = -log(rand)/S;
r(1) = 1-csi;
% b = mp/c1/(2*r(1)/(1+rmin));
t(1) = 1 - csi/(a+b*r(1));

i=1;
while t(i)>rmin
    i=i+1;
    csi = -log(rand)/S;
    r(i) = t(i-1)-csi;
%     b = mp/c1/(2*r(1)/(1+rmin));
    t(i) =  t(i-1) - csi/(a+b*r(i));

end

use = r>rmin;
r=r(use);
t=t(use);
% n=n(use);

dt = [1; t(1:end-1)] - t;
n = dt./(1+mp.*dt);


    



