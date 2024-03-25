 function   f = fun2(s,t,lambda)
 
 c = 1.5;
 c1=c+1;
 c2=c+2;
 csi = @(s,t) 1-s - (c*(t-s)./((1-t).^-c-(1-s).^-c)).^(1/(c+1));
 dcsi = @(s,t) -c^(1/c1)/c1*(t-s).^(-c/c1).*...
               ((1-t).^-c-(1-s).^-c-c*(t-s).*(1-t).^-c1)./...
               ((1-t).^-c-(1-s).^-c).^(c2/c1);
           
f = lambda*exp(-lambda*csi(s,t)).*dcsi(s,t);  
f(s==t)=lambda/2;