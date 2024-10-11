%% plot a phylogenetic tree
function [x,t,tp] = phylo(time)

s2 = 0.01;
lambda = 1;
mu = 0;
x = zeros(N,1);
t = zeros(N,1);
tp = zeros(N,1);

n = 1;
time = 0;
while n<N && n>0
    
%     p = rand;    
%     if p < lambda/(mu+lambda)
    tau = -log(rand(n,1))/lambda;
    [tb,b] = min(tau);
    dt = time+tb-t(b);
    sigma = sqrt(s2*dt);
    n = n+1; 
    x(n) = x(b) + randn*sigma;
    time = time+tb;
    t(n) = time;
    tp(n) = t(n)-t(b);
%     else 
%     tau = -log(rand(n,1))/mu;
%     [tb,b] = min(tau);
%     x(b)=[];
%     t(b)=[];
%     time = time+tb;
%     n = n-1;  
%     end

% end

end

t = time-t;


% s2 = 0.01;
% lambda = 1;
% mu = 0;
% x = zeros(N,1);
% t = zeros(N,1);
% tp = zeros(N,1);
% 
% 
% % n=0;
% % while n==0
% n = 1;
% time = 0;
% while n<N && n>0
%     
% %     p = rand;    
% %     if p < lambda/(mu+lambda)
%     tau = -log(rand(n,1))/lambda;
%     [tb,b] = min(tau);
%     sigma = sqrt(s2*tb);
%     n = n+1; 
%     x(n) = x(b) + randn*sigma;
%     time = time+tb;
%     t(n) = time;
%     tp(n) = t(n)-t(b);
% %     else 
% %     tau = -log(rand(n,1))/mu;
% %     [tb,b] = min(tau);
% %     x(b)=[];
% %     t(b)=[];
% %     time = time+tb;
% %     n = n-1;  
% %     end
% 
% % end
% 
% end
% 
% t = time-t;
    