function HypoExpPdf(lambda,zi)

%The hypoexponential is the distribution of a sum of
% k expoential random variables each with their own rate lambda


% simulaion example
% k = 7;
% lambda = 1./rand(k,1);
% R = 10000;
% z = exprnd(1,R,k)*(lambda.^-1);
% clf
% histogram(z,"Normalization","pdf");

k = length(lambda);
S = diag(-lambda);
for i=1:k-1
S(i,i+1)= lambda(i);
end
S0 = S*ones(n,1);

pdf0 = zeros(length(zi),1);
for i=1:length(zi)
Y = -expm(zi(i)*S)*S0;
pdf0(i) = Y(1);
end


hold all
plot(zi,pdf0,'linewidth',2)


