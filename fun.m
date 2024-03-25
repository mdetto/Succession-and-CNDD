function   F = fun(s,t,lambda)

c = 1.5;
 if s==t
     F = 1;
 else
     F = exp(-lambda* ...
         (1-s - (c*(t-s)./((1-t).^-c-(1-s).^-c)).^(1/(c+1))));
     F(s==t)=1;
 
 end
 
