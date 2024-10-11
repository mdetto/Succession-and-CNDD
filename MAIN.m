clear


JM = 100;
mp = [0 0.25];
c = 1.5;
c1 = c+1;
rmin = 0.5;
R = 500;
figure(1);clf
for j=1:2
csi =zeros(R*JM,1);
epsilon =zeros(R*JM,1);
N =zeros(R*JM,1);
TR =zeros(R*JM,1);
Ta =zeros(R*JM,1);
Tb =zeros(R*JM,1);
H = zeros(R,1);
n1 = 1;
for i=1:R
    i
[n,r,t,~,a,b]= LightCompetitionStrictCNDD_gost(JM,mp(j),rmin);
    
n2 = n1+length(r)-1;
csi(n1:n2,1) = a;
epsilon(n1:n2,1) = b;
N(n1:n2,1) = n;
Ta(n1:n2,1) = t;
Tb(n1:n2,1) = [1; t(1:end-1)];
TR(n1:n2,1) = r;
H(i) = length(n);
n1 = n2+1;
end

csi(n1:end,:)=[];
epsilon(n1:end,:)=[];
N(n1:end,:)=[];
TR(n1:end,:)=[];
Ta(n1:end,:)=[];
Tb(n1:end,:)=[];



rm = (1+rmin)/2;

subplot(2,2,j)
plot(csi,epsilon,'.')
% axis([1e-6 100 1e-6 100])
axis square
refline(1,0)
refline((c1-mp(j)*(1+rmin))/(c1+mp(j)*(1+rmin)),0)

% if j==1
%     mod2 = (sqrt(9-4.*csi.^2.*(c+2)./TR.*((c+2)./TR.^2.*csi-3))-3).^(1/2)...
%         ./(2*(c+2)./TR).^(1/2);
%     plot(csi,mod2,'r.')
%     pause
% end

xlabel('\xi','fontsize',15)
ylabel('\epsilon','fontsize',15)

title(['\alpha = ' num2str(mp(j))])

lambda = 2*(c+1)/(c+1+2*mp(j)*rm)*(1-rmin)/(JM+1);

Sc = (1-rmin)/lambda;

subplot(2,2,j+2)
histogram(H)
xline(Sc,'r-','linewidth',2)
xlabel('species richness')
legend('simulations','analytical solution')
axis square
pause(.1)
end
%% species richness as function of CNDD
clear

JM = [25 50 100];
mp = linspace(0,1,10);
c1 = 2.5;
rmin = 0.2;
R = 100;

S = zeros(length(mp),length(JM));
E = zeros(length(mp),length(JM));
S1 = zeros(length(mp),length(JM));
E1 = zeros(length(mp),length(JM));

for k=1:length(JM)
    for j=1:length(mp)
        j
        H = zeros(R,1);
        H1 = zeros(R,1);
        parfor i=1:R
            
            [n,r,t,~,a,b]= LightCompetitionStrictCNDD_gost(JM(k),mp(j),rmin);
            H(i,1) = length(n);
            
            % [r,n,t] = LightCompetitionStrictCNDD_linear_v2(JM(k)/(1-rmin),mp(j),rmin);
            % H1(i,1) = length(n);
        end
        
        S(j,k)=mean(H);
        E(j,k)=std(H)/sqrt(R);
        % S1(j,:)=mean(H1);
        % E1(j,1)=std(H1)/sqrt(R);
    end
    
end

figure(2);clf
for k=1:length(JM)
    % Sc = 1.2859 + (1+b*(1+rmin))*JM(k)/2;
    Sc = 1+(1+mp/c1*(1+rmin))*JM(k)/2;
    
    h = errorbar(mp,S(:,k),E(:,k),'.');hold all
    % h2 = errorbar(mp,S1,E1,'.');hold all
    plot(mp,Sc,'color',h.Color)
    xlabel('\beta','fontsize',15)
    ylabel('species richness')
    pause(.1)
end

h=get(gca,'children');
legend(h([1 3 4]),['{\itJ_M} = ' num2str(JM(3))],...
                  ['{\itJ_M} = ' num2str(JM(2))],...
                  ['{\itJ_M} = ' num2str(JM(1))])
legend('boxoff')    
set(gca,'xtick',0:0.2:1,'ytick',0:25:100)
xlim([-0.01 1.01])
axis square
%% graphical method

c = 1.5;
c1 = c+1;
d = 0.5;

% [n,r,t]= LightCompetitionStrictCNDD_gost(10,0,0.2);
   
% rp = (r.^-c - 1).^(-1/c);
% tp = (age.^-c - 1).^(-1/c);
load('GESim_v1.mat');
age = linspace(0,1.2,10000);
r1=[1.05; r(1); r(2); 0.71; r(3); r(4)];
LRS = (age./r1).^c1;

figure(1);clf
subplot(211)
plot(age,LRS([2 3 5 6],:),'k-','linewidth',2);hold all
plot(age,LRS([1 4],:),'-','linewidth',2,'color',0.8*[1 1 1])
plot(r1([2 3 5 6]),1,'ko','markersize',5,'MarkerFaceColor','k')
plot(r1([1 4]),1,'ko','markersize',5,'MarkerFaceColor',0.8*[1 1 1])
yline(1,'k--')

plot([1 1],[0 1],'r','linewidth',2)
use = age>1;
fill([age(use) fliplr(age(use))],[zeros(1,sum(use)) ones(1,sum(use))],'r','facealpha',0.1,'linestyle','none');
text(1,-0.11,'{\itt}_0')
for i=1:4
plot(r([i i]),[0 1],'r-')
plot(t([i i]),[0 (t(i)/r(i)).^c1],'r','linewidth',2)
use = age<r(i) & age>t(i);
fill([age(use) fliplr(age(use))],[zeros(1,sum(use)) fliplr(age(use)./r(i)).^c1],'r','facealpha',0.1,'linestyle','none');
text(r(i)*(1.5)^(1/c1)-0.01,1.535,num2str(i),'BackgroundColor','white')
text(r(i),-0.11,['{\itr}_' num2str(i)])
text(t(i),-0.11,['{\itt}_' num2str(i)])
end



axis([0 1.12 0 1.49])
xlabel('patch age')
ylabel('lifetime reproductive success')
set(gca,'xtick',[],'ytick',0:0.5:1)
h = get(gca,'children');
legend(h([37 35 16]),'coexisting species','excluded species','niche shadow','Location','northwest')
legend('boxoff')
%%%%%%%%%%%% dynamic simulations

J = length(r1);
imax = 81;
n = zeros(J,imax);
n(:,1)=0.1;
k = c*r1.^-c1;
t = zeros(J,1);t(1)=1;
for i=2:imax
    for j=1:J
        t(j+1) = (t(j)^-c + n(j,i-1).*k(j)).^(-1/c);
        n(j,i) = t(j)-t(j+1);
    end
end

subplot(212);
plot(0:imax-1,n([2 3 5 6],:).','-k','linewidth',2);hold all
plot(0:imax-1,n([1 4],:).','-','color',[1 1 1]*.7,'linewidth',2)
xlim([0 imax-1])
b=0;
for i=[2 3 5 6]
    b=b+1;
    text(imax,n(i,imax),num2str(b),'BackgroundColor','white')
end
set(gca,'xtick',0:20:imax-1,'ytick',0:0.1:0.3)
xlabel('simulation time step')
ylabel('population density')
sublabel

%% patch-age distribution
    
    
t = linspace(0,10,1000);

t0 =[1 2 3];
p = zeros(1000,3);
m = zeros(1000,3);
for i=1:3
p(:,i) = c/(c+1)*(1+t.^c./t0(i).^c).^((-2*c+1)/c)./t0(i);
m(:,i) = (2*c+1)*t.^(c-1)./(t0(i).^c+t.^c);
end

subplot(211)
plot(t,p)
ylabel('patch-age distribution')

legend('{\itt}_0 = 1','{\itt}_0 = 3','{\itt}_0 = 3')
legend('boxoff')

subplot(212)
plot(t,m)
xlabel('patch age')
ylabel('disturbance rate')

%% break-even time and abundance simulations
clear
rmin = 0.2;
mp = [0.0 0.75];
JM = [25 100];
rep = 100000;
c = 1.5;

R = zeros(rep*JM(2),2,2);
Z = zeros(rep*JM(2),2,2);
Z2 = zeros(rep*JM(2),2,2);
for j=1:2
    for k=1:2
      
n1 = 1;
for i=1:rep
  [j k i]

[n,r,t]= LightCompetitionStrictCNDD_gost(JM(k),mp(j),rmin);

n2 = n1+length(r)-1;
t1 = [1; t(1:end-1)];
t2 = t;
Z(n1:n2,j,k) = n.*(t1+t2)/2;
Z2(n1:n2,j,k) = c/(c-1)*(t2.^(1-c)-t1.^(1-c))./(t2.^-c-t1.^-c).*n; % exact abundance
R(n1:n2,j,k) = r;
n1 = n2+1;


end

R(n1:end,j,k)=nan;
Z(n1:end,j,k)=nan;

    end
end

save('SAD2.mat','R','Z','Z2','JM','mp','rmin')





%% break-even time distributions
clear;load('SAD2.mat')
J = 20;
MP = linspace(0,0.75,J);
  
cmap = colormap('jet');
COL = zeros(length(MP),3);
COL(:,1) = interp1(1:64,cmap(:,1),linspace(1,64,J));
COL(:,2) = interp1(1:64,cmap(:,2),linspace(1,64,J));
COL(:,3) = interp1(1:64,cmap(:,3),linspace(1,64,J));

col(1,:) =  COL(1,:);
col(2,:) =  COL(J,:);

lgnd.position{1} = [0.0635  0.5835  0.0257  0.3416];
lgnd.position{2} = [0.0586  0.1100  0.0257  0.3416];

figure(1);clf  
for k=1:2
    subplot(2,2,(k-1)*2+1)
    for i=1:J
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),MP(i),rmin);
        plot(ri,pR,'r','linewidth',2,'color',COL(i,:));hold all
    end
    xlabel('break-even time')
    caxis(MP([1 J]))
    c = colorbar;
    c.Ticks = 0:0.25:0.75;
    c.Label.String = '\alpha';
    c.Label.FontSize = 14;
    c.Label.Rotation = 0;
    c.Position = lgnd.position{k};
    pos = get(gca,'Position');
    title(['analytical results {\itJ_M} = ' num2str(JM(k))])
    for j=1:2
        
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),mp(j),rmin);
        r0i = linspace(rmin,1,50);
        xi =linspace(rmin,1-(r0i(2)-r0i(1))/2,100);
        pR0 = interp1(ri,pR(1,:),xi);
        subplot(4,2,(k-1)*4+(j-1)*2+2)
        h=histogram(R(R(:,j,k)>0,j,k),r0i,'normalization','pdf');hold all
        plot(xi,pR0,'r','linewidth',2,'color',col(j,:));hold all
        set(h,'facealpha',0.2,'facecolor',col(j,:))
        set(gca,'xtick',[0.2 1],'ytick',[0 1 2],'ylim',[0 2.5])
        t = title(['\alpha = ' num2str(mp(j),2)]);
        t.Position = [0.6000    2.1917    0.0000];
        
        if k==2 && j==1
            p = get(gca,'position');
            p(2) = pos(2)+pos(4)-p(4);
            set(gca,'position',p);
        elseif k==1 && j==2
            p = get(gca,'position');
            p(2) = pos(2);
            set(gca,'position',p);
        end
        if j==1
            set(gca,'xtick',[],'ytick',[0 1 2])
        elseif j==2
            xlabel('break-even time','fontsize',12)
        end
            
        pause(.1)
    end
end

annotation(figure(1),'textbox',...
    [0.680 0.934 0.098 0.033],'String',{'simulations'},...
    'FontSize' ,14,'linestyle','none');


%% species abundance distribution
clear;load('SAD2.mat')
J = 20;
MP = linspace(0,0.75,J);
  
cmap = colormap('jet');
COL = zeros(length(MP),3);
COL(:,1) = interp1(1:64,cmap(:,1),linspace(1,64,J));
COL(:,2) = interp1(1:64,cmap(:,2),linspace(1,64,J));
COL(:,3) = interp1(1:64,cmap(:,3),linspace(1,64,J));

col(1,:) =  COL(1,:);
col(2,:) =  COL(J,:);

lgnd.position{1} = [0.0635  0.5835  0.0257  0.3416];
lgnd.position{2} = [0.0586  0.1100  0.0257  0.3416];

figure(2);clf  
for k=1:2
    subplot(2,2,(k-1)*2+1)
    for i=1:J
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),MP(i),rmin);
        plot(Lzi,pZ(1,:),'r','linewidth',2,'color',COL(i,:));hold all
    end
    xlabel('Log abundace')
    xlim([-11.5 -0.5])
    caxis(MP([1 J]))
    c = colorbar;
    c.Ticks = 0:0.25:0.75;
    c.Label.String = '\alpha';
    c.Label.FontSize = 14;
    c.Label.Rotation = 0;
    c.Position = lgnd.position{k};
    pos = get(gca,'Position');
    title(['{\itJ_M} = ' num2str(JM(k))])
    for j=1:2
        
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),mp(j),rmin);
        Lz0i = linspace(Lzi(1),Lzi(end),50);
%         xi =linspace(rmin,1-(r0i(2)-r0i(1))/2,100);
%         pZ0 = interp1(ri,pR(1,:),xi);
        subplot(4,2,(k-1)*4+(j-1)*2+2)
        h=histogram(log(Z2(Z2(:,j,k)>0,j,k)),Lz0i,'normalization','pdf');hold all
        plot(Lzi,pZ(1,:),'r','linewidth',2,'color',col(j,:));hold all
        set(h,'facealpha',0.2,'facecolor',col(j,:))
        set(gca,'xtick',[-10 -5 0],'ytick',[0 0.2 0.4],'ylim',[0 0.4],'xlim',[-11.5 -0.5])
        t = title(['\alpha = ' num2str(mp(j),2)]);
        t.Position = [-10.0755    0.3507    0.0000];
        
        if k==2 && j==1
            p = get(gca,'position');
            p(2) = pos(2)+pos(4)-p(4);
            set(gca,'position',p);
        elseif k==1 && j==2
            p = get(gca,'position');
            p(2) = pos(2);
            set(gca,'position',p);
        end
        if j==1
            set(gca,'xtick',[],'ytick',[0 1 2])
        elseif j==2
            xlabel('Log abundance','fontsize',12)
        end
            
        pause(.1)
    end
end

annotation(figure(2),'textbox',...
    [0.680 0.956 0.098 0.033],'String',{'simulations'},...
    'FontSize' ,14,'linestyle','none');

annotation(figure(2),'textbox',...
    [0.215 0.956 0.207 0.033],...
    'String','analytical results',...
    'LineStyle','none',...
    'FontSize',14,...
    'FitBoxToText','off');

%% find intercept of S ~ CNDD

clear


S = zeros(1,100);
S0 = zeros(1,100);
JM = round(logspace(1,3,100));
R = 1000;
mp = 0.5;
rmin = .2;

for k=1:length(JM)
    
k
H = zeros(R,1);
H0 = zeros(R,1);
parfor i=1:R

% [n,r,t]= LightCompetitionStrictCNDD_gost(JM(k),mp,rmin);
% H(i,1) = length(n);
[n,r,t]= LightCompetitionStrictCNDD_linear_v2(JM(k)/(1-rmin),mp,rmin);
H0(i,1) = length(n);

end

% S(k)=mean(H);
S0(k)=mean(H0);
end


save('S_intecept2.mat','S','S0','JM','mp','rmin')

clf
subplot(121)
load S_intecept1.mat
semilogx(JM,1/2+1/2*JM.^-1);hold all
semilogx(JM,(1+mp/2.5*(1+rmin))/2+JM.^-1)
semilogx(JM,S0./JM,'bo','markersize',4)
semilogx(JM,S./JM,'ro','markersize',4)

xlabel('\itJ_M')
ylabel('s/{\itJ_M}')
legend('linearized anlytical solution','corrected solution',...
       'linearized simulations','nonlinear simulations')
title('\alpha = 0')
axis square

subplot(122)
load S_intecept2.mat
semilogx(JM,(1/2-mp/2.5)./JM + (1/2 +1/2*mp/2.5*(1+rmin)));hold all
semilogx(JM,(1+mp/2.5*(1+rmin))/2+JM.^-1)
semilogx(JM,S0./JM,'bo','markersize',4)
semilogx(JM,S./JM,'ro','markersize',4)
xlabel('\itJ_M')
% ylabel('s/{\itJ_M}')
% legend('simuations','1/2 + 1/2{\itJ_M}^-^1','1/2 + {\itJ_M}^-^1')
title('\alpha = 0.5')
axis square


%% forest structure
clear

figure(5);clf
t0 = 1;
S = [25 50 100];
mp = 0;
c = 1.5;
c1 = c+1;
rmin = 0.35;
R = 200;
g = 1;

b = 0.3;

xmax = g*1;
xmin = g*rmin^(1-b)*0.8;
x = linspace(xmin,xmax,1000);

Z1 = zeros(length(x),R);
H = zeros(R,1);
for m=1:3
for j=1:R
    j
    [n,r,t]= LightCompetitionStrictCNDD_gost(S(m),mp,rmin);
    f = g*r.^-b;
    ta = t;
    tb = [1;t(1:end-1)];
    if length(n)>1
        
        W1 = heaviside(f.*ta-x);
        W2 = heaviside(f.*tb-x);
        
        Z1(:,j) = sum(n./f.*W1 + ((x./f).^-c - tb.^-c)./(f.*c.*r.^-c1).*(1-W1).*W2);
        H(j)=length(n);
    else
        
        W1 = heaviside(r.*t-x);
        W2 = heaviside(r.*1-x);
        
        Z1(:,j) = (n.*W1 + ((x./f).^-c - tb.^-c)./(c*r.^-c1).*(1-W1).*W2)./f;
    end
    
    
end


semilogx(x,mean(Z1,2),'-');hold all
    pause(.1)
end
    
if b==0
Z0 = (1-x/g)/g;
Z0(x<g*rmin) = (1-rmin)/g;
elseif b==-1
Z0 = -log(sqrt(x/g))/g;
Z0(x<g*rmin^2) = -log(rmin)/g;
elseif b>0
Z0 = 1/(b+1)*(1-(x/g).^((b+1)/(1-b)))/g;
Z0(x<g*rmin^(1-b)) =  1/(b+1)*(1-rmin^(b+1))/g;
end

% semilogx(x,mean(Z1,2),'-');hold all
loglog(x,Z0,'-k','linewidth',2)

xlabel('tree size (x)')
ylabel('abundance')
legend(['S  = ' num2str(S(1))],...
       ['S  = ' num2str(S(2))],...
       ['S  = ' num2str(S(3))],...
        'analytical model');
    
    
%% Skewness
clear
JM = 25;
rmin = 0.2;
J = 25;
MP = linspace(0,0.75,J);
sk = zeros(J,2);
for i=1:J
    i
    [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM/(1-rmin),MP(i),rmin);
    
    % in th natural space
    zi=exp(Lzi);
    m1 = trapz(zi,pZ(1,:));
    m2 = trapz(zi,pZ(1,:)./zi.*(zi-m1).^2);
    m3 = trapz(zi,pZ(1,:)./zi.*(zi-m1).^3);
    
    sk(i,1) = m3./m2.^(3/2);   
    
    % in th log space
    m1 = trapz(Lzi,pZ(1,:).*Lzi);
    m2 = trapz(Lzi,pZ(1,:).*(Lzi-m1).^2);
    m3 = trapz(Lzi,pZ(1,:).*(Lzi-m1).^3);
    
    sk(i,2) = m3./m2.^(3/2);
end


     plot(MP,sk(:,2));hold all
     
 
 use1 = Z2(:,1,1)>0;  
 use2 = Z2(:,2,1)>0;  
sk1 = skewness(log(Z2(use1,1,1)))
sk2 = skewness(log(Z2(use2,2,1)))
abs(1-sk2/sk1)*100

use1 = Z2(:,1,2)>0;  
use2 = Z2(:,2,2)>0;  
sk1 = skewness(log(Z2(use1,1,2)))
sk2 = skewness(log(Z2(use2,2,2)))
abs(1-sk2/sk1)*100  

 
 use1 = Z2(:,1,1)>0;  
 use2 = Z2(:,2,1)>0;  
sk1 = skewness(Z2(use1,1,1))
sk2 = skewness(Z2(use2,2,1))
abs(1-sk2/sk1)*100

use1 = Z2(:,1,2)>0;  
use2 = Z2(:,2,2)>0;  
sk1 = skewness(Z2(use1,1,2))
sk2 = skewness(Z2(use2,2,2))
abs(1-sk2/sk1)*100  

%% adatptive evolution
R = 100;
J = 55;
Tmax = round(2.^linspace(4,8,J));
mp = [0 0.25 0.5];
p = 0.05;
s = ones(R,J,length(mp));


for k=1:length(mp)
    TR0 = cell(R,1);
    parfor i=1:R
        [k i]
        s(i,:,k) = AdpaptEvol_v3(mp(k),p,Tmax);
    end
end

clf
col(1,:) = [0         0.4471    0.7412];
col(2,:) = [0.8510    0.3255    0.0980];
col(3,:) = [0.4667    0.6745    0.1882];
qt = 0.45;
Lti = linspace(4,7.5,100);
for i=1:3

t0 = repmat(Tmax,R,1); 
s0 = s(:,:,i);
x = (s0(:));
y = log2(t0(:));

pdf0 = ksdensity(y,Lti);
pdf1 = ksdensity(y,Lti,'weights',x);
pdf2 = ksdensity(y,Lti,'weights',x.^2);

mi = pdf1./pdf0*mean(x);
vi = pdf2./pdf0*mean(x.^2) - mi.^2;

subplot(211)
gray_area(Lti,mi,mi-sqrt(vi/R),mi+1*sqrt(vi/R),'color',col(i,:))
lgnd{4-i} = ['\alpha = ' num2str(mp(i))];
end


subplot(211)
set(gca,'xtick',4:8,'xticklabel',2.^(4:8),...
        'ytick',0:10:50,'xlim',[4 7.5])
h = get(gca,'children');
legend(h([1 3 5]),lgnd)
legend('boxoff')
xlabel('evolutionary time (speciation events)')
ylabel('species richness')

Tstar = (1-rmin)/2/(0.1*sqrt(2/pi))/(0.01*2);

%%  adatptive evolution distributions

rmin = 0.2;
R = 200;
J = 75;
Tmax = logspace(1,3.8,J);
mp = [0 0.1 0.2];
S0 = [20 40];
TR1 = zeros(R,S0(1),length(mp));
TR2 = zeros(R,S0(2),length(mp));

for k=1:length(mp)
parfor i=1:R
    [k i]
    
    [r,rm,T0] = AdpaptEvol_v2(rmin,mp(k),0,S0(1),'Smax');
    TR1(i,:,k)=r(1:S0(1));
    
    [r,rm,T0] = AdpaptEvol_v2(rm,mp(k),T0,S0(2),'Smax');
    TR2(i,:,k)=r(1:S0(2));

end

end

%% plot analysis above
col(1,:) = [0         0.4471    0.7412];
col(2,:) = [0.8510    0.3255    0.0980];
col(3,:) = [0.4667    0.6745    0.1882];

subplot(223);cla
ri = linspace(rmin+.01,1-.01,100);
N = zeros(length(ri),3);
for i=1:3
    TR0 = TR1(:,:,i);
N(:,i) = ksdensity(TR0(:),ri,'support',[rmin 1]);
end

for j=1:3
plot(ri,N(:,j),'color',col(j,:),'linewidth',2);hold all
end
title(['{\its} = ' num2str(S0(1))])
ylim([0 12])
xlabel('break-even time')
ylabel('pdf')

subplot(224);cla
ri = linspace(rmin+.01,1-.01,50);
N = zeros(length(ri),3);
for i=1:3
    TR0 = TR2(:,:,i);
N(:,i) = ksdensity(TR0(:),ri,'support',[rmin 1]);
end
for j=1:3
plot(ri,N(:,j),'color',col(j,:),'linewidth',2);hold all
end
title(['{\its} = ' num2str(S0(2))])
ylim([0 12])
xlabel('break-even time')
 
%% adpative evolution - species richness
clear 
R = 100;
Tmax = 150;
mp = linspace(0,0.5,7);

s0 = ones(length(mp),R);
sm = ones(length(mp),R);

h0 = ones(length(mp),R);
% hm = ones(length(mp),R);

p = 0.05;
sigma = 0.03;

for k=1:length(mp)
    parfor i=1:R
        [k i]
        [s0(k,i),sm(k,i)] = AdpaptEvol_v3(mp(k),p,0,Tmax);
%         [h0(k,i),hm(k,i)] = AdpaptEvol_v4(mp(k),p,3*Tmax);

    end
end
 
subplot(212)
% plot(mp,mean(sm,2),'ro','markersize',4,'MarkerFaceColor','r');hold all
plot(mp,mean(s0,2),'ro','markersize',4,'MarkerFaceColor','r');hold all

plot([mp; mp],mean(s0,2).'+[std(s0,0,2) -std(s0,0,2)].'/sqrt(R),'r-')


% plot(mp,mean(h0,2),'o')
lsline
xlabel('\beta')
ylabel('species richness')
legend('single trait evolution','multi trait evolution')
legend('metacommunity','local community')
legend('boxoff')

%% single vs.multi-trait selection

R=100;SM = 20:10:70;K=length(SM);
s0 = ones(R,K);
sm = ones(R,K);
ts = ones(R,K);
h0 = ones(R,K);
hm = ones(R,K);
th = ones(R,K);

p = 0.05;
Tmax = 20000;
for k=1:6
    k
    parfor i=1:R
        [s0(i,k),sm(i,k),~,~,~,ts(i,k)] = AdpaptEvol_v3(0,p,SM(k),Tmax);
        [h0(i,k),hm(i,k),~,~,~,th(i,k)] = AdpaptEvol_v4(0,p,SM(k),Tmax);
    end
end

clf
plot(SM+0.1,mean(s0),'ro','markersize',4,'MarkerFaceColor','r');hold all
plot(SM-0.1,mean(h0),'bo','markersize',4,'MarkerFaceColor','b')
lsline


% plot([SM; SM]+0.1,mean(s0)+ 2*[-1 1]'*std(s0)/sqrt(R),'b-')
% plot([SM; SM]-0.1,mean(h0)+ 2*[-1 1]'*std(h0)/sqrt(R),'r-')

plot([SM; SM]+0.1,quantile(s0,[0.25 0.75]),'r-')
plot([SM; SM]-0.1,quantile(h0,[0.25 0.75]),'b-')

legend('single trait evolution','multi trait evolution')
% legend('metacommunity','local community')
legend('boxoff')

xlabel('metacommunity size ({\itJ_M})')
ylabel('species richness')

set(gca,'xlim',[18 72],'ylim',[0 60],'xtick',0:20:70,'ytick',0:20:70)
axis square

%% abundance as function of break-even time (simulations)
clear
JM = 100;
beta = [0 0.4 0.8];
rmin = 0.2;
imax = 1000;
ri = linspace(rmin,1,101);
M = zeros(length(ri)-1,3);
for j=1:3

N =zeros(imax*JM,2);
R =zeros(imax*JM,2);
n1 = 1;
for i=1:imax
    i
[n,r,t,~,a,b]= LightCompetitionStrictCNDD_gost(JM,beta(j),rmin);
    
n2 = n1+length(r)-1;
t1 = [1; t(1:end-1)];
t2 = t;

R(n1:n2) = r;
% N(n1:n2) = c/(c-1)*(t2.^(1-c)-t1.^(1-c))./(t2.^-c-t1.^-c).*n; % exact abundance
N(n1:n2)=n;
n1 = n2+1;
end


N(n1:end)=[];
R(n1:end)=[];

for w=1:length(ri)-1
    use = R>ri(w)&R<=ri(w+1);
    M(w,j)=mean(N(use));
end
% pdf0(:,1) = ksdensity(R,r,'Support',[rmin 1]);
% pdf1(:,1) = ksdensity(R,r,'weights',N,'Support',[rmin 1]);
% pdf2(:,1) = ksdensity(R,r,'weights',N.^2,'Support',[rmin 1]);
% 
% mi(:,j) = pdf1./pdf0*mean(N);
% vi(:,j) = pdf2./pdf0*mean(N.^2) - mi(:,j).^2;
end

% col(1,:) = [0         0.4471    0.7412];
% col(2,:) = [0.8510    0.3255    0.0980];
% col(3,:) = [0.4667    0.6745    0.1882];

% figure(1);clf
% for j=1:3
% gray_area(y,mi(:,j),mi(:,j)-sqrt(vi(:,j)/imax),mi(:,j)+sqrt(vi(:,j)/imax),'color',col(j,:))
% hold all
% end
% plot(y,mi);hold all
% xlabel('break-even time','fontsize',15)
% ylabel('abundance','fontsize',15)
% pause(.1)

plot(ri(1:100)+(ri(2)-ri(1))/2,M);
xlabel('break-event time')
ylabel('initial recruits')


%% abundance as function of break-even time (analytical)
JM = 100;
beta = [0 0.3 0.8];
rmin = 0.2;
c=1.5;
c1=2.5;
lambda = JM/(1-rmin);

[M0,ti,R_r,r,pZ,Lzi,s,pN] = RenewalProcess2(lambda,0,rmin,1.5);
R = zeros(length(r),length(beta));
R(:,1) = R_r;
for j=2:length(beta)
[M0,ti,R_r,r,pZ,Lzi,s,pN] = RenewalProcess2(lambda,beta(j),rmin,1.5);
 R(:,j) = R_r;
end
rp = (r.^-c - 1).^(-1/c);

% r = linspace(rmin,1,1001);
% rp = (r.^-c - 1).^(-1/c);
% R = zeros(length(r),length(beta));
% 
% R(:,1) = 2*(1-exp(-lambda*(1-r)))/lambda;
% for j=1:length(beta)
%  h = lambda.*(1/2/beta(j)+r/c1);
%  k = lambda*(1-r);
%  R(:,j) = (1-exp(-k) - h.*exp(h).*...
%      (expint(h) - expint(k+h)))/beta(j);
% end
subplot(211)
plot(r,R)
xlabel('break-event time')
ylabel('initial recruits')
legend([repmat('\beta = ',3,1) num2str(beta.')])
legend('Boxoff')

subplot(212)
plot(rp,R)
xlabel('break-event time')
ylabel('initial recruits')
legend([repmat('\beta = ',3,1) num2str(beta.')])
legend('Boxoff')

h(1).XData=r;
h(2).XData=r;
h(3).XData=r;

h(1).YData=R(:,3);
h(2).YData=R(:,2);
h(3).YData=R(:,1);

%% sensitivity to c and rmin
clear;clf
rmin = 0.2;
beta = 0.5;
JM = 100;
c = linspace(1,2,20);
J = length(c);
cmap = colormap('jet');ncmap = length(cmap);
COL = zeros(length(c),3);
COL(:,1) = interp1(1:ncmap,cmap(:,1),linspace(1,ncmap,J));
COL(:,2) = interp1(1:ncmap,cmap(:,2),linspace(1,ncmap,J));
COL(:,3) = interp1(1:ncmap,cmap(:,3),linspace(1,ncmap,J));

subplot(221)
for j=1:J

     [M,ti,pR,ri,pZ,Lzi,SL(j),pN] = RenewalProcess(JM/(1-rmin),beta,rmin,c(j));
        plot(Lzi,pZ(1,:),'r','linewidth',2,'color',COL(j,:));hold all

end
xlabel('distribution')
ylabel('log abundance')
h=colorbar;
set(h,'Ticks',0:0.5:1,'TickLabels',[1 1.5 2],'Position',[0.1641    0.8071    0.0135    0.0851])
set(h.Label,'string','\itc')

set(h,'ytick')
subplot(223)
plot(c,SL,'linewidth',1)
xlabel('parameter c')
ylabel('species richness')

rmin = linspace(0.2,0.8,J);
subplot(222)
for j=1:J

     [M,ti,pR,ri,pZ,Lzi,SL(j),pN] = RenewalProcess(JM/(1-rmin(j)),beta,rmin(j),c);
        plot(Lzi,pZ(1,:),'r','linewidth',2,'color',COL(j,:));hold all

end

ylabel('log abundance')
colorbar
h=colorbar;
set(h,'Ticks',0:0.5:1,'TickLabels',[0.2 0.5 0.8],'Position',[0.6073    0.8071    0.0135    0.0851])
set(h.Label,'string','{\itr}_{min}')

subplot(224)
plot(rmin,SL,'linewidth',1)
xlabel('parameter r_{min}')

%% adptive evolution SAD
