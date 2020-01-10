function [prob_class] = script_3
% Imprecise likelihoods, (standard) Bayesian approach
%
% B. Ristic, RMIT University, March 2018

% compute the confusion matrix
confusion_matrix_prob  = [1 1/2 1/3; 0 1/2 1/3; 0 0 1/3] ;


feature_vec = [1 1 1 1 1 1 2 1 1 1 1 2 3 1 1 1 1];
L = length(feature_vec);
prior = [1/3 1/3 1/3];
prob_class = Bayesian(confusion_matrix_prob,prior,feature_vec);
%     % Run Bayesian possibilistic classifier

%disp(prob_class(1,:))
h=figure(1);
set(h,'Position',[50 50 400 300]);
plot([0:L],prob_class(1,:),'r','Linewidth',2);
axis([0 L -0.1 1.1]); 
hold on
plot([0:L],prob_class(2,:),'b-o');
plot([0:L],prob_class(3,:),'g--','Linewidth',2);
hold off;
legend('Class 1','Class 2','Class 3');
ylabel('Probability of a class');
xlabel('Measurement index k');
%title('Test3: Bayesian');



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


