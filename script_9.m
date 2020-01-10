function post_matrix = script_9(diagnosticity, n, prior)
% Testing scenario 1 - TBM solution (posterior tables) 
% n - the size of the confusion matrix
% diagnosticity: 
% prior_m - prior probability in form of a bba m (as a vector)
% Example: script_9(5, 3, [2/5 2/5 1/5]');

addpath 'TBM'
addpath 'TBM\FMT'

x = diagnosticity/(diagnosticity + n - 1);
element_off_diagonal = (1-x)/(n-1);

% mapping from 1:n indices to 1:2^n indices of bba
seq = gen_seq(n,2);
cnt = 0;
for i=1:2^n
    if sum(seq(i,:))==1
        cnt = cnt + 1;
        mapp(cnt) = i; 
    end
end
% Confusion matrix (each colum 2^n elements: bba on feature space)

confusion_matrix = zeros(2^n,n);
prior_m = zeros(2^n,1);
for i=1:n
    confusion_matrix(mapp(i),:) = ones(1,n)*element_off_diagonal;
    confusion_matrix(mapp(i),i) = x; 
    prior_m(mapp(i)) = prior(i);
end
confusion_matrix
prior_m

% % modification: least commited bbas
% for i=1:n
%     mm = pign2m([confusion_matrix([2 3 5],i)]);
%     confusion_matrix(:,i) = mm;
% end


% Compute the posterior matrix
for j=1:n
    [mH] = GBT1(confusion_matrix,mapp(j));
     m = conjun(mH,prior_m);
     p = pignistic(m);
     post_matrix(j,:) = p(mapp)';
end
post_matrix'

end
%%
