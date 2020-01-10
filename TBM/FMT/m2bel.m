function [Bel] = m2Bel(m)
% computing FMT from m to bel
% in = m vector
% out = bel vector

bel = b2bel(m2b(m));
Bel = bel / (1-m(1));