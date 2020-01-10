function script_15(num_mc)
% Application of Imprecise Probability Theory to solve
% a modified Testing Scenario 4.
% The confusion matrix is given by intervals.
%
% Input:
% num_MC - number of Monte Carlo runs
%
% Output: 
%   IPclass2: probability interval for class 2
%   ProbClass2: precise probability (Bayesian, correct, for comparison)
% 


% define the confusion matrix parameters
d1=5;
n=3;
prior = [2/5 2/5 1/5];

% compute the confusion matrix for measurement generation
diag_elem1 = d1/(d1 + n - 1);
element_off_diagonal1 = (1-diag_elem1)/(n-1);
confusion_matrix1 = element_off_diagonal1*ones(n,n);
for i=1:n
    confusion_matrix1(i,i) = diag_elem1;
end
%confusion_matrix1

class = 2; % Set the true class to be 2
M=20;   % total number of measurements (duration of scenario)

if isfile('res.mat')
    load res.mat;
    start = i;
    fprintf('starting from i=%d   ... press "Enter" to continue\n',i);
    pause;
else
    start=0;
end

for i=start+1:num_mc
    % generate input data (vector of measurements)
    feature_vec = NaN* ones(M,1);
    for m = 1:M
        % use confusion matrix 1 to generate data
        feature_vec(m) = resample(confusion_matrix1(:,class),1);
    end
    % Run Bayesian probabilistic classifier (correct & mismatched)
    res = Bayesian(confusion_matrix1,prior,feature_vec);
    ProbClass2All(i,:) = res(2,:);
    
    % Run IP classifier (imprecise confuion matrix)
    [Results] = IPclass(feature_vec, prior,1);  %
    UPclass2All_a(i,:) = Results(4,:);
    LPclass2All_a(i,:) = Results(3,:);
 
    save res i ProbClass2All LPclass2All_a UPclass2All_a ;
    
    ProbClass2 = mean(ProbClass2All,1);
    
    LPclass2a = mean(LPclass2All_a,1);
    UPclass2a = mean(UPclass2All_a,1);
    %IPclass2a = [LPclass2a; UPclass2a];
    
    if i > 1
        figure(50);
        kk = [0:M];
        plot(kk,ProbClass2,'b-',kk,LPclass2a,'r-.',kk,UPclass2a,'r--');
        axis([0 M 0.  1.02]); 
        hold off;
        xlabel('Measurement index k');
        ylabel('Probability of class 2');
        legend('Bayesian','Lower Prob','Upper Prob');
        drawnow;
    end
end
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
function feature_vec = gen_meas(M,true_class)
% M - number of features

feature_vec = zeros(1,M);
d1 = 5;
n=3;
diag_elem1 = d1/(d1 + n - 1);
element_off_diagonal1 = (1-diag_elem1)/(n-1);
confusion_matrix1 = element_off_diagonal1*ones(n,n);
for i=1:n
    confusion_matrix1(i,i) = diag_elem1;
end

for m = 1:M
    % use confusion matrix 1 to generate data
    feature_vec(m) = resample(confusion_matrix1(:,true_class),1);
end

end

%%
%-----------------------------------------------------------------
% The  classifier code
function [Results] = IPclass(feature_vec, prior,type)
% Imprecise likelihoods, imprecise probability
% Testing scenario 4 (modified)
% Input: feature vector, prior pmf, and 'type' determines the 
%        confusion matrix (1 or 2)
% A. Benavoli, modified by B Ristic.


%%
p0=prior;
P=[p0]; %initial credal set for class x_i

%% Conditional credal sets Z|x given as inequality constraints.

%variables of linear programming problem
%g(z_1|x_1),g(z_2|x_1),....g(z_2|x_3),g(z_3|x_3) 9 variables

alpha = [1 1 1; -1 -1 -1];
beta = [1; -1];
Aineq = [kron(eye(3),alpha); kron(eye(9),beta)];

if type == 1
    gamma = [0.75; -0.65];
    delta = [0.20; -0.10];
else
    gamma = [0.75; -0.7];
    delta = [0.17; -0.12];
end
bineq=[beta' beta' beta' gamma' delta' delta' delta' gamma' delta' delta' delta' gamma']';

lower_bnd=zeros(9,1);

%% Inference Loop
h = waitbar(0,'Running IP code');
Results=zeros(3*2,length(feature_vec));
for obs=1:length(feature_vec) %loop on observations
    waitbar(obs/length(feature_vec),h);
    observedZ=feature_vec(obs); %observation at this time instant
    
    LP=zeros(3,size(P,1)); %this vector will include lower probability of x_i
    UP=zeros(3,size(P,1)); %this vector will include upper probability of x_i
    
    for ip=1:size(P,1) %run over vertices credal set K(X)
        p=P(ip,:);
        for class_x=1:3 %loop on class x_i
            w=zeros(1,3);
            
            %computing lower probability of classes
            w(class_x)=1;
            myfun = @(nu) abs(solve_LP(p,Aineq,bineq,lower_bnd,observedZ,w,nu)); %this solves the min poblem
            optnu=fminbnd(myfun,0,1);%this implements bisection on nu
            LP(class_x,ip)=optnu;
                
           
            %computing upper probability of classes
            w=ones(1,3)-w;
            myfun = @(nu) abs(solve_LP(p,Aineq,bineq,lower_bnd,observedZ,w,nu)); %this solves the min poblem
            optnu=fminbnd(myfun,0,1);%this implements bisection on nu
            UP(class_x,ip)=1-optnu;
            
        end
        
    end
    LP=min(LP,[],2);
    UP=max(UP,[],2);
    Res=[LP(1) UP(1) LP(2) UP(2) LP(3) UP(3)]';
    Res=round(Res,3);
    Results(:,obs)=Res;
    
    
    %Computing posterior credal set of p(x_i) after this observation
    %This is actually an outer approximation of the credal set but it is
    %exact in this example
    Acc=[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
    bcc=[Res(2,1);-Res(1,1);Res(4,1);-Res(3,1);Res(6,1);-Res(5,1)];
    Aeq=[1 1 1];
    beq=1;
    Vert=[];
    for i=1:100
        c=randn(1,3);%generating random search direction
        sol=linprog(c,Acc,bcc,Aeq,beq,zeros(3,1));
        Vert=[Vert;sol'];
    end
    
    P=unique(Vert,'rows');
    disp(Results)

    
end
Results=[[p0(1) p0(1) p0(2) p0(2) p0(3) p0(3)]', Results];
close(h);

if 1 < 0
% plotting
L = length(feature_vec);
figure(20);
plot([0:L],Results(1,:),'r-','Linewidth',2);
hold on
axis([0 L -0.1 1.1]);
plot([0:L],Results(2,:),'r--','Linewidth',2);
legend('LP','UP')
ylabel('Lower/Upper Probability of class 1');
xlabel('Measurement index k');
hold off

figure(21);
plot([0:L],Results(3,:)+0.003,'b-','Linewidth',2); %+0.003 is only for plotting reasons
hold on
axis([0 L -0.1 1.1]);
plot([0:L],Results(4,:)+0.003,'b-.','Linewidth',2);
legend('LP','UP')
ylabel('Lower/Upper Probability of class 2');
xlabel('Measurement index k');
hold off

figure(22);
plot([0:L],Results(5,:),'g-','Linewidth',1);
hold on
axis([0 L -0.1 1.1]);
plot([0:L],Results(6,:),'g:','Linewidth',2);
legend('LP','UP')
ylabel('Lower/Upper Probability of class 3');
xlabel('Measurement index k');
hold off;
end


end

%%
function funval=solve_LP(p,A,b,lower_bnd,observedZ,w,nu)
%constraint: credal sets

Z=observedZ;


if Z==1
    c=[(fun(w,1)-nu)*p(1) 0 0  (fun(w,2)-nu)*p(2) 0 0 (fun(w,3)-nu)*p(3) 0 0];
    
elseif Z==2
    c=[ 0 (fun(w,1)-nu)*p(1) 0  0 (fun(w,2)-nu)*p(2) 0  0 (fun(w,3)-nu)*p(3) 0];
    
elseif Z==3
    c=[ 0 0 (fun(w,1)-nu)*p(1) 0 0 (fun(w,2)-nu)*p(2) 0  0 (fun(w,3)-nu)*p(3)];
end


[sol,funval]=linprog(c,A,b);

end

%% function of interest
function val = fun(w,arg)
val=w(1)*Indicator(1,arg)+w(2)*Indicator(2,arg)+w(3)*Indicator(3,arg);
end

%%
% indicator function
function val = Indicator(v,cond)
val=0;
if v==cond
    val=1;
end
end




