% Script_13: Testing scenario 5
% B. Ristic, copyright RMIT University, 2018

addpath 'TBM'
addpath 'TBM\FMT'

% TBM approach
m1 = [0 0.6 0.4 0]';  % bba (probability) on {x1,x2}
m1e = balloon_ext(m1,0,1);
m2 = [0 0.7 0.3 0]';  % bba (probability) on {x2,x3}
m2e = balloon_ext(m2,1,0);

m12e = conjun(m1e,m2e);
[prob1] = pignistic(m12e);
[prob2] = voorbraak(m12e);  % Table 11 results
prob2

% Probabilistic approach (beta = 1)
p1u = [0.6 0.4 1];
p1 = p1u/sum(p1u);   % probability function 1
p2u = [1 0.7 0.3];
p2 = p2u/sum(p2u);   % probability function 2

p = p1 .* p2 / sum(p1 .* p2);   % same as Voorbraak result!!!
p
% possibilistic approach (beta = 1)
beta = 1;
pi1 = [0.6 0.4 beta]/max([0.6 0.4 beta]);
pi2 = [beta 0.7 0.3]/max([beta 0.7 0.3]);

pi_temp = pi1 .* pi2;
pi = pi_temp / max(pi_temp);
prob = pi / sum(pi);   % same as Voorbraak result!!!
prob


