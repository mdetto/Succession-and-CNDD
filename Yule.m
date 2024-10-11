% simulate Yule process with Ornstein–Uhlenbeck process
function [r,t] = Yule(tmax,lambda,mu,sigma,theta,x0)



% theta = 1;
% mu = 0;
% sigma = 0.3;
% lambda = 1;
% x0 = 10;
S = 100;
c = 1.5;
mu1 = log((mu.^-c-1).^(-1/c));
x01 = log((x0.^-c-1).^(-1/c));

x = zeros(S,1);
t = zeros(S,1);

x(1) = x01;
n = 1;
time = 0;
while time<tmax
    
    
    tau = -log(rand(n,1))/lambda;
    [tb,b] = min(tau);
    time = time+tb;
    n = n+1;
        
        % BM
        % sigma = sqrt(s2*(time-t(b)));
        % x(n) = x(b) + randn*sigma;

        % O-U
        T = time-t(b);
        K2 = 1/2*sigma^2/theta*(1-exp(-2*theta*T));
        x(n) = exp(-theta*T)*x(b) + exp(-theta*T)*(exp(theta*T)-1)*mu1 + sqrt(K2)*randn;     
        t(n) = time;
    
end

x = x(2:n);
t = t(n) - t(2:n);
r = (exp(x).^-c+1).^(-1/c);



