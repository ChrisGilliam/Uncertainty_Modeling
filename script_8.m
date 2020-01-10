function [prob_c,prob_m,poss2p_m] = script_8(d1, n, prior, num_mc)
% Testing scenario with incorrect confsuion matrix: 
% Comparing via Monte Carlo the Bayesian vs Possibilistic
% Input:
% diagnosticity - diagonal / off-diagonal element of Conf. matrix
% n - the size of the confusion matrix
% prior- prior probability vector
% num_MC - number of Monte Carlo runs
%
% Example:  script_8(5, 3, [2/5 2/5 1/5],1000);
% B. Ristic, RMIT University, March 2018

% compute the confusion matrix for measurement generation
diag_elem1 = d1/(d1 + n - 1);
element_off_diagonal1 = (1-diag_elem1)/(n-1);
confusion_matrix1 = element_off_diagonal1*ones(n,n);
for i=1:n
    confusion_matrix1(i,i) = diag_elem1;
end
confusion_matrix1

confusion_matrix2 = [0.42 0.375 0.15; 0.18 0.4 0.3; ...
                    1-(0.42+0.18) 1-(0.375+0.4) 1-(0.15+0.3)];
confusion_matrix2

class = 2; % Set the true class to be 2
M=25;   % total number of measurements (duration of scenario)
for i=1:num_mc
    % generate input data (vector of measurements)
    feature_vec = NaN* ones(M,1);
    for m = 1:M
        % use confusion matrix 1 to generate data
        feature_vec(m) = resample(confusion_matrix1(:,class),1);
    end
    % Run Bayesian probabilistic classifier (correct & mismatched)
    prob_class_c(:,:,i) = Bayesian(confusion_matrix1,prior,feature_vec);
    prob_class_m(:,:,i) = Bayesian(confusion_matrix2,prior,feature_vec);
    % Run Bayesian possibilistic classifier (correct & mismatched)
    [pi,pi2p]  = Jeremie(confusion_matrix1,prior,feature_vec);
    poss_class_c(:,:,i) = pi;
    poss2p_class_c(:,:,i) = pi2p;
    [piA,pi2pA] = Jeremie(confusion_matrix2,prior,feature_vec);
    poss_class_m(:,:,i) = piA;
    poss2p_class_m(:,:,i) = pi2pA;
end
prob_c = mean(prob_class_c,3);
prob_m = mean(prob_class_m,3);
poss_c = mean(poss_class_c,3);
poss_m = mean(poss_class_m,3);
poss2p_c = mean(poss2p_class_c,3);
poss2p_m = mean(poss2p_class_m,3);

figure(20);
plot([0:M],prob_c(2,:),'bv-',[0:M],prob_m(2,:),'gs:',...
     [0:M],poss2p_m(2,:),'r--');
axis([0 M 0.35 1.05]); 
hold off;
xlabel('Measurement index k');
ylabel('Probability of class 2');

legend('Correct model','Model-mismatch, Bayesian ',...
    'Model-mismatch, Possibilistic', 'Location','East');

title('testing scenario 2');

end
%%
function resample_idx= resample(w,L)
% resampling 
% function resample_idx= resample(w,L)
% w- the weights with sum(w)= 1
% L- no. of samples you want to resample
% resample_idx- indices for the resampled particles

resample_idx= [];
[notused,sort_idx]= sort(-w);   %sort in descending order
rv= rand(L,1);
i= 0; 
threshold= 0;
while ~isempty(rv),
    i= i+1; 
    threshold= threshold+ w(sort_idx(i));
    rv_len= length(rv);
    idx= find(rv>threshold); 
    resample_idx= [ resample_idx; sort_idx(i)*ones(rv_len-length(idx),1) ];
    rv= rv(idx);
end;
end


%%
function class_prob = Bayesian(confusion_matrix,prior,feature_vec)
M = length(feature_vec);
n = size(confusion_matrix,1);
class_prob(:,1) = prior';
for m=1:M    
    for j=1:n
        class_prob(j,m+1) = class_prob(j,m)*confusion_matrix(feature_vec(m),j);
    end
    norm_const = sum(class_prob(:,m+1));
    class_prob(:,m+1) = class_prob(:,m+1)/norm_const;
end

end

%%
function [poss_class,poss2p_class] = Jeremie(confusion_matrix,prior,feature_vec)
M = length(feature_vec);
n = size(confusion_matrix,1);
% Convert probabilities to possibilities
poss_prior = prior/max(prior);
poss_conf_mat = nan*ones(n,n);
for i=1:n
    poss_conf_mat(:,i) = confusion_matrix(:,i)/max(confusion_matrix(:,i));
end
%
% Class possibilities
poss_class(:,1) = poss_prior';
poss2p_class(:,1) = poss_prior' / sum(poss_prior);
for m=1:M    
    for j=1:n
        poss_class(j,m+1) = poss_class(j,m)*poss_conf_mat(feature_vec(m),j);
    end
    norm_const = max(poss_class(:,m+1));
    poss_class(:,m+1) = poss_class(:,m+1)/norm_const;
    %
    poss2p_class(:,m+1) = poss_class(:,m+1) / sum(poss_class(:,m+1));
end

end
