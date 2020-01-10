function [Pl] = m2Pl(m)
% computing FMT from m to pl
% in = m vector
% out = pl vector

pl = q2pl(m2q(m));
Pl = pl / (1-m(1));
