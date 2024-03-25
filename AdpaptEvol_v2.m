function [r,rm,T] = AdpaptEvol_v2(r0,alfa,T,Tmax,type)

opts = optimset('TolX',1e-18,'TolFun',1e-18);

p = 0.5;
rmin = 0.2;
c = 1.5;
c1 = c+1;
lambda = 0.1;
sigma = .1;
r=r0;
rm = r0;
if strcmp(type,'Tmax')
while T<Tmax
    Nm = length(rm);
    dt = -log(rand)/lambda;
    j=randi(Nm);
    %         dr = sigma*randn*sqrt(dt);
    dr = sigma*randn;
    rnew = rm(j) + dr;
    if rnew>rmin
    rm = [rm; rnew];
    r   = zeros(length(rm),1);
    t   = 1;
    i=0;
    while min(rm)<t
        i=i+1;
        r(i) = max(rm(rm<t));
        t = fminbnd(@fun2,0,r(i),opts);
    end
    r = r(r>rmin);

    % randomly keep excluded
    member = ismember(rm,r);
    np = sum(~member);
    member(~member) = rand(np,1)>p;
    rm=rm(member);
    end
    T=T+dt;
end

%% stop at Smax
elseif strcmp(type,'Smax')
    N=1;
    while N<Tmax
    Nm = length(rm);
    j=randi(Nm);
    %         dr = sigma*randn*sqrt(dt);
    dr = sigma*randn;
    rnew = rm(j) + dr;
    if rnew>rmin
    rm = [rm; rnew];
    r   = zeros(length(rm),1);
    t   = 1;
    i=0;
    while min(rm)<t
        i=i+1;
        r(i) = max(rm(rm<t));
        t = fminbnd(@fun2,0,r(i),opts);
    end
    r = r(r>rmin);

    % randomly keep excluded
    member = ismember(rm,r);
    np = sum(~member);
    member(~member) = rand(np,1)>p;
    rm=rm(member);
    end
    N=length(r);
    end
end
%% function for t
    function y = fun2(ti)

        ni = (ti^-c - t^-c)*r(i)^c1/c;
        Ri = ni/(1-alfa*ni);
        y = abs(Ri - t + ti);
    end

end



