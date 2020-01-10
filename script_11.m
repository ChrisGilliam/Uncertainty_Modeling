function [pign_prob_class] = script_11
% Testing scenario 3: TBM
%
% B. Ristic, RMIT University, April 2018



addpath 'TBM'
addpath 'TBM\FMT'

prior_m = zeros(8,1);
prior_m(end) = 1;

conf_matrix_bba = zeros(8,3);
conf_matrix_bba(2,1) = 1;
conf_matrix_bba(4,2) = 1;
conf_matrix_bba(8,3) = 1;

mapp = mapping_fun(3);

feature_vec = [1 1 1 1 1 1 2 1 1 1 1 2 3 1 1 1 1];
L = length(feature_vec);

pign_prob_class = tbm_class(conf_matrix_bba,prior_m,mapp,feature_vec);

h=figure(1);
set(h,'Position',[50 50 400 300]);
plot([0:L],pign_prob_class(1,:),'r','Linewidth',2);
axis([0 L -0.1 1.1]); 
hold on
plot([0:L],pign_prob_class(2,:),'b-o');
plot([0:L],pign_prob_class(3,:),'g--','Linewidth',2);
hold off;
legend('Class 1','Class 2','Class 3','Location','northwest');
ylabel('Probability of a class');
xlabel('Measurement index k');
title('Test3: TBM');



% 
% % Plotting
% figure(22);
% errorbar([1:M+1],prob(1,:),prob_err(1,:));
% axis([0 M -0.1 1.1]);
% figure(23);
% errorbar([1:M+1],prob(2,:),prob_err(2,:));
% axis([0 M -0.1 1.1]);
% figure(24);
% errorbar([1:M+1],prob(3,:),prob_err(3,:));
% axis([0 M -0.1 1.1]);
% legend('show');

end

%%
function class_prob = tbm_class(confusion_matrix,prior_m,mapp,feature_vec)
M = length(feature_vec);
%n = length(mapp);
mH_all = prior_m;
p = pignistic(mH_all);
class_prob(:,1) = p(mapp)';
for m=1:M  
    [mH] = GBT1(confusion_matrix,mapp(feature_vec(m)));
     mH_all = conjun(mH,mH_all);
     p = pignistic(mH_all);
     class_prob(:,m+1) = p(mapp)';
     % TBM has a problem unless the bba is normalised
     mH_all = normalise(mH_all);
end

end


%%
function mapp = mapping_fun(n)
% mapping from 1:n indices to 1:2^n indices of bba
seq = gen_seq(n,2);
cnt = 0;
for i=1:2^n
    if sum(seq(i,:))==1
        cnt = cnt + 1;
        mapp(cnt) = i; 
    end
end
end