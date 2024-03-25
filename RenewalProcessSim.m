%% Simulate a renewal process

function t = RenewalProcessSim(a,b,tmax)

t  = zeros(1000,1);

i=1;
while t(i)<tmax
    i=i+1;
    t(i) =  t(i-1) -log(rand)/(a+b*t(i-1));

end

t([1 i:end])=[];



    



