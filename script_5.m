function post_poss_matrix = script_5(diagnosticity, n, prior)
% Possibilistic classifier, test scenario 1
%
% Input: n - the size of the confusion matrix
%        diagnosticity
%        prior_probability
% Example: scrip_3(5, 3, [2/5 2/5 1/5]);


x = diagnosticity/(diagnosticity + n - 1);

element_off_diagonal = (1-x)/(n-1);
confusion_matrix = element_off_diagonal*ones(n,n);
for i=1:n
    confusion_matrix(i,i) = x;
end

poss_prior = prior/max(prior);

poss_conf_mat = nan*ones(n,n);
for i=1:n
    poss_conf_mat(:,i) = confusion_matrix(:,i)/max(confusion_matrix(:,i));
end
poss_conf_mat   % Table 2.(a)

for i=1:n
    for j=1:n
        post_poss_matrix(j,i) = poss_conf_mat(i,j)*prior(j);
    end
    post_poss_matrix(:,i) = post_poss_matrix(:,i) /  max(post_poss_matrix(:,i));
end
post_poss_matrix  % Table 2.(b)

% Convert possibilities to probabilities (normalise)
for i=1:n
    post_poss2prob_mat(:,i) = post_poss_matrix(:,i)/sum(post_poss_matrix(:,i));
end
post_poss2prob_mat   % Table 1.(b)



