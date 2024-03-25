%% Simulation of the Ornstein–Uhlenbeck process
function x = OU_sim(x0,t,sig,th,mu)

N = length(t);

W = cumsum(sqrt(diff(exp(2*th*t))).*randn(1,N-1));
ex = exp(-th*t(2:N));

x(1) = x0;
x(2:N) = x0*ex+mu*(1-ex)+sig*ex.*W/sqrt(2*th);

