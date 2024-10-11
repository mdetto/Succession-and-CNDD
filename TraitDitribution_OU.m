function P = TraitDitribution_OU(ti,tmax,lambda,mu,sigma,theta,x0,type)

c = 1.5;
mu1 = log((mu.^-c-1).^(-1/c));
x01 = log((x0.^-c-1).^(-1/c));

M0 = @(t) x01*exp(-theta*t) + exp(-theta*t).*(exp(theta*t)-1)*mu1;
V0 = @(t) 1/2*sigma^2/theta*(1-exp(-2*theta*t));
P = zeros(1,length(ti));

xi = -1/c*log(ti.^-c-1);

if strcmp(type,'pdf')

for i=1:length(xi)
    if ti(i)==0 || ti(i)==1
    P(1,i) = 0;
    else
    P(1,i) = integral(@(t) lambda*exp(-lambda*t)./(1-exp(-lambda*tmax)) ... 
        .*1./sqrt(2*pi*V0(tmax-t)).*exp(-(xi(i)-M0(tmax-t)).^2./(2*V0(tmax-t)))...
        .*ti(i).^-(c+1)./(ti(i).^-c-1),0,tmax);
    end
end

elseif strcmp(type,'cdf')
for i=1:length(xi)
    if ti(i)==0
    P(1,i) = 0;
    elseif ti(i)==1
       P(1,i) = 1;
    else
    P(1,i) = integral(@(t) lambda*exp(-lambda*t)./(1-exp(-lambda*tmax)) ... 
        .*1/2.*erfc(-(xi(i)-M0(tmax-t))./sqrt(2*V0(tmax-t))),0,tmax);
    end
end

else
    disp('error: type not defined properly')
end


