function  script_10(diagnosticity, n, prior,num_mc)
%  Testing scenario 2 - TBM solution (recursive, Monte Carlo runs) 
% diagnosticity - diagonal / off-diagonal element of Conf. matrix
% n - the size of the confusion matrix
% prior- prior probability vector
% num_MC - number of Monte Carlo runs
%
% >> script_10(5, 3,[2/5 2/5 1/5]',100);
% B. Ristic, RMIT University, April 2018
% 

addpath 'TBM'
addpath 'TBM\FMT'

diag_elem = diagnosticity/(diagnosticity + n - 1);
element_off_diagonal = (1-diag_elem)/(n-1);
prob_confusion_matrix = element_off_diagonal*ones(n,n);
for i=1:n
    prob_confusion_matrix(i,i) = diag_elem;
end
prob_confusion_matrix
% mapping from 1:n indices to 1:2^n indices of bba
seq = gen_seq(n,2);
cnt = 0;
for i=1:2^n
    if sum(seq(i,:))==1
        cnt = cnt + 1;
        mapp(cnt) = i; 
    end
end
% Confsuion matrix (each colum 2^n elements: bba on feature space)

confusion_matrix = zeros(2^n,n);
prior_m = zeros(2^n,1);
for i=1:n
    confusion_matrix(mapp(i),:) = ones(1,n)*element_off_diagonal;
    confusion_matrix(mapp(i),i) = diag_elem; 
    prior_m(mapp(i)) = prior(i);
end

confusion_matrix
prior_m
% % modification: least commited bbas
% for i=1:n
%     mm = pign2m([confusion_matrix([2 3 5],i)]);
%     confusion_matrix(:,i) = mm;
% end


class = 2; % Set the true class to be 2
M=25;   % total number of measurements (duration of scenario)
for i=1:num_mc
    % generate input data (vector of measurements)
    feature_vec = NaN* ones(M,1);
    for m = 1:M
        feature_vec(m) = resample(prob_confusion_matrix(:,class),1);
    end
    % Run TBM classifier (output is pignistic class probability)
    pign_prob_class(:,:,i) = tbm_class(confusion_matrix,prior_m,mapp,feature_vec);
    %m0(i,:) = mass_emptyset(confusion_matrix,prior_m,mapp,feature_vec);
    % For comparison, we also run Bayesian probabilistic classifier
    prob_class(:,:,i) = Bayesian(prob_confusion_matrix,prior,feature_vec);
end

% PLOTTING
mean_pign_prob = mean(pign_prob_class,3);
pign_prob_err = std(pign_prob_class,0,3);

h=figure(20);
set(h,'Position',[450 50 400 300]);
errorbar([1:M+1],mean_pign_prob(1,:),pign_prob_err(1,:),'-s','Markersize',6, ...
    'MarkerFaceColor','red');
axis([0 M -0.1 1.1]);
hold on;
errorbar([1:M+1],mean_pign_prob(2,:),pign_prob_err(2,:),'-s','Markersize',6, ...
    'MarkerFaceColor','green');
hold off;


% report the difference between Bayesian and TBM
h=figure(21);
set(h,'Position',[650 350 400 300]);
mean_prob = mean(prob_class,3);
plot([1:M+1],abs(mean_pign_prob(2,:)-mean_prob(2,:)) ,'r','Linewidth',2);

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
function class_prob = tbm_class(confusion_matrix,prior_m,mapp,feature_vec)
M = length(feature_vec);
%n = length(mapp);
mH_all = prior_m;
class_prob(:,1) = prior_m(mapp)';
for m=1:M  
    [mH] = GBT1(confusion_matrix,mapp(feature_vec(m)));
     mH_all = conjun(mH,mH_all);
     p = pignistic(mH_all);
     class_prob(:,m+1) = p(mapp)';
     % Normalise bbm if the mass at zero too high
     if mH_all(1) > 0.9
        mH_all = normalise(mH_all);
     end
end

end

%%
function mass0 = mass_emptyset(conf_matrix_tbm,prior_tbm,mapp,feature_vec)

M = length(feature_vec);
mass0 = NaN*ones(1,M+1);
mass0(1) = 0;
mH_all = prior_tbm;
for m=1:M  
    [mH] = GBT1(conf_matrix_tbm,mapp(feature_vec(m)));
     mH_all = conjun(mH,mH_all);
    mass0(m+1) = mH_all(1);
end

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