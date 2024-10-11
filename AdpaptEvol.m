function [r,t,T] = AdpaptEvol(r0,t,alfa,T,Tmax,type)


opts = optimset('TolX',1e-18,'TolFun',1e-18);

rmin = 0.2;
c = 1.5;
c1 = c+1;
lambda = 0.01;
sigma = .1;
r = r0;
N = length(r);
% for s=1:K

if strcmp(type,'Tmax')
    k = 0; C=0;
while T<Tmax

        tau = -log(rand(N,1))/lambda;
        [dt,j]=min(tau);
%         dr = sigma*randn*sqrt(dt);
        dr = sigma*randn;
        rnew = r(j) + dr;

       if rnew>rmin && rnew<1
           [r,t] = LightCompetition(r,t,rnew);
           
%            if ismember(rnew,r)
%            k=k+1;
%            C = C + sum(dr>0);
%            end
       end
       N = length(r);
       
       T=T+dt;
       
%     plot(r,1,'ko')
%     xlim([rmin 1])
%     title(['s = ' num2str(N)])
%     pause(.1)
end       

% p = C/k;

elseif strcmp(type,'Smax')
    
    while N<Tmax

        tau = -log(rand(N,1))/lambda;
        [dt,j]=min(tau);
        rnew = r(j) + sigma*randn*sqrt(dt);
       if rnew>rmin && rnew<1
           [r,t] = LightCompetition(r,t,rnew);
       end
       N = length(r);
       
       T=T+dt;
       
%     plot(r,1,'ko')
%     xlim([rmin 1])
%     title(['s = ' num2str(N)])
%     pause(.1)
    end   
end
        
%% select species 
function [r0,t] = LightCompetition(tr,t,rnew)
        

    tr = [tr; rnew];
    r0  = zeros(length(tr),1);
    
    if length(tr)>2 && max(tr)>rnew
    i = sum(tr>rnew);
    r0(1:i) = tr(tr>rnew);
    
    else
        
    r0(1) = max(tr(tr<1));
    t(1) = fminbnd(@fun1,0,r0(1),opts);
    i=1;
    end
    
    while min(tr)<t(i)
        i=i+1;

        r0(i) = max(tr(tr<t(i-1)));
        t(i)  = fminbnd(@fun2,0,r0(i),opts);

    end
    use = r0>rmin;
    r0=r0(use);
    t=t(use);
  

%% function for i=1
function y = fun1(t1)

    n1 = (t1^-c - 1)*r0(1)^c1/c;
    R1 = n1/(1-alfa*n1);
    y = abs(R1 - 1 + t1);

end

%% function for i>1
function y = fun2(ti)
       
    ni = (ti^-c - t(i-1)^-c)*r0(i)^c1/c;
    Ri = ni/(1-alfa*ni);
    y = abs(Ri - t(i-1) + ti);
end

end

end


