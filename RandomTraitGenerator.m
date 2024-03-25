function tr = RandomTraitGenerator(N,tmax,lambda,mu,sigma,theta,x0)

c = 1.5;
mu1 = log((mu.^-c-1).^(-1/c));
x01 = log((x0.^-c-1).^(-1/c));

M0 = @(t) x01*exp(-theta*t) + exp(-theta*t).*(exp(theta*t)-1)*mu1;
V0 = @(t) 1/2*sigma^2/theta*(1-exp(-2*theta*t));


z =  rand(N,1);
xr = zeros(N,1);
for i=1:N
    N-i
    xr(i) = fzero(@fun,2);
end

c = 1.5;
tr = (exp(xr).^-c+1).^(-1/c);


function y = fun(xi)

y = integral(@(t) lambda*exp(-lambda*t)./(1-exp(-lambda*tmax)) ... 
        .*1/2.*erfc(-(xi-M0(tmax-t))./sqrt(2*V0(tmax-t))),0,tmax) - z(i);
end
end