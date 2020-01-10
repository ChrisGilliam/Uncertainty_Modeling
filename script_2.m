function  prob_class = script_2(diagnosticity, n, prior, num_mc)
% Bayesian classifier, testing scenario 1
% Monte Carlo results (true class is x2)
% Input:
% diagnosticity - diagonal / off-diagonal element of Conf. matrix
% n - the size of the confusion matrix
% prior- prior probability vector
% num_MC - number of Monte Carlo runs
%
% Example:  script_2(5, 3, [2/5 2/5 1/5],100);
% B. Ristic, RMIT University, March 2018

% compute the confusion matrix
diag_elem = diagnosticity/(diagnosticity + n - 1);
element_off_diagonal = (1-diag_elem)/(n-1);
confusion_matrix = element_off_diagonal*ones(n,n);
for i=1:n
    confusion_matrix(i,i) = diag_elem;
end
confusion_matrix


class = 2; % Set the true class to be 2
M=25;   % total number of measurements (duration of scenario)
for i=1:num_mc
    % generate input data (vector of measurements)
    feature_vec = NaN* ones(M,1);
    for m = 1:M
        feature_vec(m) = resample(confusion_matrix(:,class),1);
    end
    % Run Bayesian probabilistic classifier
    prob_class(:,:,i) = Bayesian(confusion_matrix,prior,feature_vec);
end
prob = mean(prob_class,3);
prob_err = std(prob_class,0,3);

ind1 = find(prob_err(1,:)+prob(1,:)> 1);
pos_bar(1,:) = prob_err(1,:);
pos_bar(1,ind1) = 1-prob(1,ind1);
ind2 = find(prob_err(2,:)+prob(2,:)> 1);
pos_bar(2,:) = prob_err(2,:);
pos_bar(2,ind2) = 1-prob(2,ind2);

ind3 = find(-prob_err(1,:)+prob(1,:)< 0);
neg_bar(1,:) = prob_err(1,:);
neg_bar(1,ind3) = prob(1,ind3);
ind4 = find(-prob_err(2,:)+prob(2,:)< 0);
neg_bar(2,:) = prob_err(2,:);
neg_bar(2,ind4) = prob(2,ind4);

% Plotting
figure(22);
%errorbar([0:M],prob(1,:),prob_err(1,:),'-s','Markersize',6, ...
%    'MarkerFaceColor','red');
errorbar([0:M],prob(1,:),neg_bar(1,:),pos_bar(1,:),'-s','Color','red','Markersize',6, ...
    'MarkerFaceColor','red');
xlabel('Measurement index k');
ylabel('Probability of a class');
%axis([1 M+1 -0.1 1.1]);
hold on
errorbar([0:M],prob(2,:),neg_bar(2,:),pos_bar(2,:),'-s','Color','blue','Markersize',6, ...
    'MarkerFaceColor','blue');
%axis([1 M+1 -0.1 1.1]);
hold off
legend('Class 1','Class 2', 'Location','East');
axis([0 M  -0.05 1.05]);

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
function poss_class = Jeremie(confusion_matrix,prior,feature_vec)
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
for m=1:M    
    for j=1:n
        poss_class(j,m+1) = poss_class(j,m)*poss_conf_mat(feature_vec(m),j);
    end
    norm_const = max(poss_class(:,m+1));
    poss_class(:,m+1) = poss_class(:,m+1)/norm_const;
end
end
