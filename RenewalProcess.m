%% Compute the break-even time distribution and species abundance distribution
% in the local community given a uniformely distribuited
% metacommunity. Specifically, inter-arrival times are exponential
% random variables with rate S. 

function [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(S,mp,rmin)

bin = 100;
tmax = 1-rmin;

if mp==0

%% renwal function
ti = linspace(0,2*tmax,250);
M = S/2*ti;
use = ti>tmax;
M(use) = (S*tmax + 1 - exp(-S*(ti(use)-tmax)))/2;


%% species trait distribution
ri = linspace(rmin,1,bin);
pR = S*(exp(-S*(1-ri))+1)/(S*tmax + 1 - exp(-S*tmax));
  

%% SAD
zmax = 0.5;
Lzi = linspace(-12,log(zmax),bin);
zi = exp(Lzi);

SL  = (S*tmax + 1 - exp(-S*tmax))/2;

use=zi> 2*tmax*(1-tmax);
xup = zi./(1-tmax);
xup(use)=1+sqrt(1-2*zi(use));
xlo=1-sqrt(1-2*zi);

pZ = S/2*zi.*exp(-S/2*(1-sqrt(1-2*zi)))./sqrt(1-2*zi)/SL + ...
      S^2/(4*SL)*(expint(S*xlo/2)-expint(S*xup/2)).*zi;

pZ(end) = 0;


%% SAD 2

tup = tmax - zi/(2*(1-tmax));
zcr = 2*(1-tmax)^2;
tup(zi>zcr) =  1 - sqrt(2*zi(zi>zcr));

x1    = @(t,z) 1-t-(sqrt((1-t).^2-2*z));
x2    = @(t,z) 1-t+(sqrt((1-t).^2-2*z));
dxdz = @(t,z) 1./(sqrt((1-t).^2-2*z));

f    = @(x) S/2.*exp(-S/2.*x);

pi2 = zeros(1,length(zi));
for i=1:length(zi)-1
    
       pi2(i) = integral(@(t) ...
         zi(i).*f(x1(t,zi(i))).*dxdz(t,zi(i)).*S/2,0,tup(i));

end


p1 = zi.*f(1-sqrt(1-2*zi))./(sqrt(1-2*zi));
p1(end)=0;


pZ(2,:) = (p1+pi2)/SL;

%% P(N)

pN = zi./SL.*f(zi).*(1+S*(tmax-zi/2)/2);

%% E(Z|t)
ti = linspace(rmin,1,1000);
EZt = zeros(1,length(ti));
for i=1:length(ti)
    
       EZt(i) = integral(@(x) S*exp(-S*x).*(2*x.*(ti(i)-x)),0,ti(i)-rmin);
           
end
%  plot(ti,EZt)

%%%%%%%   CASE WITH CNDD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else  

c1 = 2.5;
tmax = 1-rmin;
a = 1/2 + mp/c1;
b = mp/c1;
% c0 = 1/(integral(@(t) exp(-S*a./(1+b*t).*t).*(a-b*t),0,tmax));
% SL = (2*S^2*(S-a*c0)+b*c0*(1-2*S))/(2*S^3).*(1-exp(-S*tmax)) + c0*(a-1/2*b*tmax)*tmax;

ap = a-b/(2*S);
SL = (1-ap + b*(1-2*S)/(2*S^2)).*(1-exp(-S*tmax)) + S*(ap-1/2*b*tmax)*tmax;

%% renwal function
ti = linspace(0,tmax/(a-b*tmax),250);

M = S*(ap-b/2*ti).*ti;



%% species trait distribution
ri = linspace(rmin,1,bin);

% pR = c0/SL*(S/c0 - a - b/S*(1-1/2/S))*exp(-S*(1-ri)) + c0/SL*(a-b*(1-ri));
pR = S/SL*(1 - ap - b/S*(1-1/2/S))*exp(-S*(1-ri)) + S/SL*(ap-b*(1-ri));

%% SAD



zmax = (1+mp-sqrt(1+2*mp))/mp^2;
Lzi = linspace(-12,log(zmax),bin);
zi = exp(Lzi);

zcr = tmax/2*(2*a-tmax*(1+2*b))/(a-b*tmax)/(a-(b-mp)*tmax);
use = zi<zcr;
xup = 1-mp*zi+sqrt((1-mp*zi).^2-2*zi);
xup(use) = (1 - tmax - mp*zi(use) - sqrt((1 - tmax - mp*zi(use)).^2 + 4*b*(1-tmax)*zi(use))) ...
           ./(2*b*(tmax-1));
xdn = 1-mp*zi-sqrt((1-mp*zi).^2-2*zi);

t = @(x,z) 1-mp*z-z./x-x/2;
f = @(t,x) S*(a-b*t)./(1+b*x).^2.*exp(-S*(a-b*t)./(1+b*x).*x);
% m = @(t) c0*(a-b*t);
m = @(t) S*(ap-b*t);

pi = zeros(1,length(zi));
for i=1:length(zi)-1
    
    pi(i) = integral(@(x) zi(i).*f(t(x,zi(i)),x).*m(t(x,zi(i))).*(mp*x+1)./x,xdn(i),xup(i));
end

p0 = zi.*f(0,1-mp*zi-sqrt((1-mp*zi).^2-2*zi)).*...
    ((mp*(1-mp*zi)+1)./(sqrt((1-mp*zi).^2)-2*zi) - mp);
p0(end)=0;

% this part is negligible
% p02 = zi.*f(0,1-mp*zi+sqrt((1-mp*zi).^2-2*zi)).*...
%     ((mp*(1-mp*zi)+1)./(sqrt((1-mp*zi).^2)-2*zi) + mp);
% p02(zi<zcr)=0;
% p0 = p0+p02;

pZ = (p0+pi)/SL;
pZ = pZ/trapz(Lzi,pZ);
%% SAD 2 - integrate in t
k = a-b*tmax;
tup = tmax - ...
    k*(tmax-1+mp*zi + sqrt((tmax-1+mp*zi).^2+2*zi*(2*k-1)))./(2*k-1);

zcr = ((a-b*tmax-1)+sqrt((1-a+b*tmax)^2+2*mp*(1-tmax))).^2/(2*mp^2);
tup(zi>zcr) = 1 -   sqrt(2*zi(zi>zcr)) - mp*zi(zi>zcr);

x    = @(t,z) 1-t-mp*z-(sqrt((1-t-mp*z).^2-2*z));
dxdz = @(t,z) (mp*(1-t-mp*z)+1)./sqrt((1-t-mp*z).^2-2*z)-mp;
f    = @(t,x) S*(a-b*t)./(1+b*x).^2.*exp(-S*(a-b*t)./(1+b*x).*x);

pi = zeros(1,length(zi));
for i=1:length(zi)-1
    
       pi(i) = integral(@(t) zi(i).*...
           f(t,x(t,zi(i))).*dxdz(t,zi(i)).*m(t),0,tup(i));
           
end


x    = @(t,z) 1-t-mp*z+(sqrt((1-t-mp*z).^2-2*z));
dxdz = @(t,z) (mp*(1-t-mp*z)+1)./sqrt((1-t-mp*z).^2-2*z)+mp;
use = find(zi>zcr);
for i=use(1:end-1)
    
       pi(i) = pi(i) + integral(@(t) zi(i).*...
           f(t,x(t,zi(i))).*dxdz(t,zi(i)).*m(t),0,tup(i));
end


pZ(2,:) = (p0+pi)/SL;




%% P(N)
pN = zeros(size(zi));
for i=1:length(zi)
pN(i) = zi(i)./SL.*(f(0,zi(i)) + integral(@(t) f(t,zi(i)).*m(t),0,tmax-k*zi(i)));
end

% g = 1+b*zi;
% www = (-exp(-a*S.*zi./b).*(g.^2 + (g+a.*S.*zi).^2) + ...
%   exp(-S.*zi.*k).*...
%   (g.^2 + (g+S*k*zi.*g).^2))./...
%   (b*g.*S^2.*zi.^3); 
% 
% pN(2,:) = zi./SL.*(f(0,zi) + c0*www);

end

