function [out] = q2pl(in)
% computing FMT from q to pl.
% in = q vector
% out = pl vector
in(1) = 0;
out = abs(b2m(in));
out = out';
