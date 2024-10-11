function [xi,yi,pdf0,pdf1] = KernelRegression(x,y,xi)

%% Nadaraya–Watson kernel regression
use = ~isnan(y) & ~isnan(x);
if nargin<3
qt = quantile(x(use),[0.025 0.975]);
xi = linspace(qt(1),qt(2),25);
end
pdf0 = ksdensity(x(use),xi);
pdf1 = ksdensity(x(use),xi,'weights',y(use))*mean(y(use));

yi = pdf1./pdf0;
if nargout==0
    hold on
    plot(xi,yi)
end
