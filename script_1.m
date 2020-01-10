function script_1(diagnosticity, n, prior)
% Bayesian classifier, test scenario 1
% Input: n - the size of the confusion matrix
%        diagnosticity
%        prior_probability
% Example: matlab_scrip_1(5, 3, [2/5 2/5 1/5]);


x = diagnosticity/(diagnosticity + n - 1);

element_off_diagonal = (1-x)/(n-1);
confusion_matrix = element_off_diagonal*ones(n,n);
for i=1:n
    confusion_matrix(i,i) = x;
end
confusion_matrix

for i=1:n
    denom = sum(confusion_matrix(i,:).*prior);
    for j=1:n
        posterior_matrix(j,i) = confusion_matrix(i,j)*prior(j)/denom;
    end
end
posterior_matrix