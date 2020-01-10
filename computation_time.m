function elapsed_time = computation_time
elapsed_time = zeros(1,5);
% Bayesian approach (standard)
tic;
script_3;
script_3;script_3;script_3;script_3;
elapsed_time(1)=toc/5;

% Mahler-Bayesian
tic;
script_4;script_4;script_4;script_4;script_4;
elapsed_time(2)=toc/5;

% Possibilistic
tic;
script_7;script_7;script_7;script_7;script_7;
elapsed_time(3)=toc/5;

% TBM
tic;
script_11;script_11;script_11;script_11;script_11;
elapsed_time(4)=toc/5;

% Imprecise prob
tic;
script_14;
elapsed_time(5)=toc;


