function m2 = normalise(m1);
%
m2 = m1/(1-m1(1));
m2(1) = 0;