function pk = PK(ki,S,kmax)
% Compute the probability distribution of K in the local community
% for a uniform dustribution of K in the metacommunity
% the distribution does not include species 1
% P_ks1(k) is the distribution of kstar1
% inputs: 
%           ki: a vecor of K when the distribution will be computed
%           S : species richnes in the metacommunity
%           Kmax: maximum valu of K

c  = 1.5;
k0 = c+1;
lambda = (S+1)/kmax;

pk = zeros(size(ki));
for i=1:length(ki)

if ki(i)==kmax
    dk = (kmax-k0)/1000;
    pk(i) = integral(@(k) 1./(kmax - k).*P_ks1(k),k0,ki(i)-dk) +dk*lambda/(S+1);
else
    pk(i) = integral(@(k) 1./(kmax - k).*P_ks1(k),k0,ki(i));
end
end

function y = P_ks1(k)
    
t = ((c+1)./k).^(1/(c+1));
L = (c+1).*t.^-c./(c-c*t+1) - k0;
dL_dt = c*(c+1)./(t.^c.*(c-c*t+1)).*(1./(c-c*t+1)-1./t);
dt_dk = -(c+1)^(-c/(c+1))./k.^((c+2)/(c+1)); 
y = dL_dt.*dt_dk.*lambda.*exp(-lambda.*L);

end

end