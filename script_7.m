function [poss_class_prob] = script_7
% Imprecise likelihood - Possibilistic classifier
%
% B. Ristic, RMIT University, March 2018

% compute the confusion matrix
confusion_matrix_poss  = [1 1 1; 0 1 1; 0 0 1];

feature_vec = [1 1 1 1 1 1 2 1 1 1 1 2 3 1 1 1 1];
L = length(feature_vec);
poss_prior = [1 1 1];

poss_class = Jeremie(confusion_matrix_poss,poss_prior,feature_vec);
poss_class_prob = poss_class ./sum(poss_class);


h=figure(1);
set(h,'Position',[50 50 400 300]);
plot([0:L],poss_class(1,:),'r','Linewidth',2);
axis([0 L -0.1 1.1]); 
hold on
plot([0:L],poss_class(2,:),'b-o');
plot([0:L],poss_class(3,:),'g--','Linewidth',2);
hold off;
legend('Class 1','Class 2','Class 3','Location','North');
ylabel('Possibility of a class');
xlabel('Measurement index k');
%title('Test 3: Possibilistic');

if 1 < 0
h=figure(2);
set(h,'Position',[450 50 400 300]);
plot([0:L],poss_class_prob(1,:),'r','Linewidth',2);
axis([0 L -0.1 1.1]); 
hold on
plot([0:L],poss_class_prob(2,:),'b-o');
plot([0:L],poss_class_prob(3,:),'g--','Linewidth',2);
hold off;
legend('Class 1','Class 2','Class 3','Location','North');
ylabel('Probability of a class');
xlabel('Measurement index k');
%title('Test 3: Possibilistic');
end


end


%%
function poss_class = Jeremie(poss_conf_mat,poss_prior,feature_vec)
M = length(feature_vec);
n = size(poss_conf_mat,1);

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
