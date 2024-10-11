%% Simulate birth-death (extinction) process with Browniwn motion
function [r,t] = Yule_V2(tmax,lambda,mu,sigma,theta,x0,S)



% theta = 1;
% mu = 0;
% sigma = 0.3;
% lambda = 1;
% x0 = 10;
c = 1.5;
mu1 = log((mu.^-c-1).^(-1/c));
x01 = log((x0.^-c-1).^(-1/c));


t = zeros(S,1);


s0 = x01/10;
x(:,1) = randn(S,1)*s0+x01;
n = 1;
time = 0;
while time<tmax
    

    tau = -log(rand(S,1))/lambda;
    [tb,b] = min(tau);
    j = randi(S);
    time = time+tb;
%     n = n+1;
        
        % BM
        % sigma = sqrt(s2*(time-t(b)));
        % x(n) = x(b) + randn*sigma;

        % O-U
        T = time-t(b);
        K2 = 1/2*sigma^2/theta*(1-exp(-2*theta*T));
        x(j) = exp(-theta*T)*x(b) + exp(-theta*T)*(exp(theta*T)-1)*mu1 + sqrt(K2)*randn;     
        t(j) = time;
    
end

r = (exp(x).^-c+1).^(-1/c);