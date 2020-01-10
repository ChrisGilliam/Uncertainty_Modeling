function m = conjun(m1,m2)
% Conjunctive rule of combination

if length(m1) ~= length(m2), error('cannot combine');end

two_n = length(m1);
n = log2(two_n);

M = 1;
for i=1:n
    M = [M zeros(2^(i-1)); M M];
end
J = zeros(two_n);
for i=1:two_n
    J(i,two_n-i+1)=1;
end
JM = J * M;
JMJ = JM * J;
q1 = JMJ * m1;
q2 = JMJ * m2;
q = q1 .* q2;
m = J * inv(M) * J * q;


   