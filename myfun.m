function y = myfun(mp,ks1,kmax,S0)

c = 1.5;
c1 = c+1;
y = zeros(size(ks1));
for i=1:length(ks1)
    
y(i) = integral(@(K) (2./(2*mp*c*c1^(-(2*c+1)/c1)*K.^(-1/c1)+1).*kmax/(S0+1)).^-1,ks1(i),kmax);

end