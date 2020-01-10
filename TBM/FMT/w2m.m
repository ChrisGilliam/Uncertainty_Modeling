function [out] = w2m(in)
% computing FMT from w to m. 
% in = w vector
% out = m vector
out = q2m(w2q(in))';
