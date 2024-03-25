function P = OU_pdf(ti,tmax,lambda,mu,sigma,theta,x0)


M0 = @(t) x0*exp(-theta*t) + exp(-theta*t).*(exp(theta*t)-1)*mu;
V0 = @(t) 1/2*sigma^2/theta*(1-exp(-2*theta*t));
P = zeros(1,length(ti));


for i=1:length(ti)
    P(1,i) = integral(@(t) lambda*exp(-lambda*t)./(1-exp(-lambda*tmax)) ... 
        .*1./sqrt(2*pi*V0(tmax-t)).*exp(-(ti(i)-M0(tmax-t)).^2./(2*V0(tmax-t))),0,tmax);
    end
end