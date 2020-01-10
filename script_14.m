function [prob_class] = script_14
% Imprecise likelihoods, imprecise probability
%
% A. Benavoli, IDSIA, Sept. 2018

feature_vec = [1 1 1 1 1 1 2 1 1 1 1 2 3 1 1 1 1];

%%
p0=[1/3 1/3 1/3]; 
P=[p0]; %initial credal set for class x_i

%% Conditional credal sets Z|x given as inequality constraints.

%variables of linear programming problem
%g(z_1|x_1),g(z_2|x_1),....g(z_2|x_3),g(z_3|x_3) 9 variables

Aineq=[1 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 1 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 1 1 1
    -1 0 0 0 0 0 0 0 0;
    0 -1 0 0 0 0 0 0 0;
    0 0 -1 0 0 0 0 0 0;
    0 0 0 -1 -1 0 0 0 0;
    0 0 0 0 0 -1 0 0 0;
    0 0 0 0 0 0 -1 -1 -1];
bineq=[1;0;0;1;0;1;-1;0;0;-1;0;-1];
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
    for i=1:50
        c=randn(1,3);%generating random search direction
        sol=linprog(c,Acc,bcc,Aeq,beq,zeros(3,1));
        Vert=[Vert;sol'];
    end
    
    P=unique(Vert,'rows');
    disp(Results)

    
end
Results=[[p0(1) p0(1) p0(2) p0(2) p0(3) p0(3)]', Results];
close(h);

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

%%
function funval=solve_LP(p,A,b,lower_bnd,observedZ,w,nu);
%constraint: credal sets

Z=observedZ;


if Z==1
    c=[(fun(w,1)-nu)*p(1) 0 0  (fun(w,2)-nu)*p(2) 0 0 (fun(w,3)-nu)*p(3) 0 0];
    
elseif Z==2
    c=[ 0 (fun(w,1)-nu)*p(1) 0  0 (fun(w,2)-nu)*p(2) 0  0 (fun(w,3)-nu)*p(3) 0];
    
elseif Z==3
    c=[ 0 0 (fun(w,1)-nu)*p(1) 0 0 (fun(w,1)-nu)*p(2) 0  0 (fun(w,3)-nu)*p(3)];
end


[sol,funval]=linprog(c,A,b,[],[],lower_bnd);

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


