function [poss] = script_6(d, n, prior, num_mc)
% Possibilistic classifier via Monte Carlo
% Input:
% diagnosticity - diagonal / off-diagonal element of Conf. matrix
% n - the size of the confusion matrix
% prior- prior probability vector
% num_MC - number of Monte Carlo runs
%
% Example:  script_6(5, 3, [2/5 2/5 1/5],5000);
% B. Ristic, RMIT University, March 2018

% compute the confusion matrix for measurement generation
diag_elem = d/(d + n - 1);
element_off_diagonal = (1-diag_elem)/(n-1);
confusion_matrix = element_off_diagonal*ones(n,n);
for i=1:n
    confusion_matrix(i,i) = diag_elem;
end
confusion_matrix
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class = 2; % Set the true class to be 2
M=25;   % total number of measurements (duration of scenario)
for i=1:num_mc
    % generate input data (vector of measurements)
    feature_vec = NaN* ones(M,1);
    for m = 1:M
        feature_vec(m) = resample(confusion_matrix(:,class),1);
    end
    % Run Bayesian possibilistic classifier
    [pi,pi2p]  = Jeremie(confusion_matrix,prior,feature_vec);
    poss_class(:,:,i) = pi;
    poss2p_class(:,:,i) = pi2p;
end
poss = mean(poss_class,3);
poss_err = std(poss_class,0,3);

ind1 = find(poss_err(1,:)+poss(1,:)> 1);
pos_bar(1,:) = poss_err(1,:);  % positive bar
pos_bar(1,ind1) = 1-poss(1,ind1);
ind2 = find(poss_err(2,:)+poss(2,:)> 1);
pos_bar(2,:) = poss_err(2,:);
pos_bar(2,ind2) = 1-poss(2,ind2);

ind3 = find(-poss_err(1,:)+poss(1,:)< 0);
neg_bar(1,:) = poss_err(1,:);   % negative bar
neg_bar(1,ind3) = poss(1,ind3);
ind4 = find(-poss_err(2,:)+poss(2,:)< 0);
neg_bar(2,:) = poss_err(2,:);
neg_bar(2,ind4) = poss(2,ind4);

h1=figure(21);
set(h1,'Position',[250 50 400 300]);
errorbar([0:M],poss(1,:),neg_bar(1,:),pos_bar(1,:),'-s','Color','red','Markersize',6, ...
    'MarkerFaceColor','red');
%axis([0 M -0.1 1.1]);
hold on; 
errorbar([0:M],poss(2,:),neg_bar(2,:),pos_bar(2,:),'-s','Color','blue','Markersize',6, ...
    'MarkerFaceColor','blue');
hold off;
legend('Class 1','Class 2', 'Location','East');
axis([0 M  -0.05 1.05]);

poss2p = mean(poss2p_class,3);
%poss_err = std(poss_class,0,3);
h2=figure(22);
set(h2,'Position',[450 50 400 300]);
%plot([1:M+1],prob(1,:),'rs-','Markersize',8);
hold on;
plot([1:M+1],poss2p(1,:),'rs-','Linewidth',2);
%plot([1:M+1],prob(2,:),'gs-','Markersize',8)
plot([1:M+1],poss2p(2,:),'bs-','Linewidth',2);
axis([0 M -0.05 1.05]);
hold off;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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
