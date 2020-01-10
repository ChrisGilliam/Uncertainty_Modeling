function [out] = q2b (in)
% computing FMT from q to b. Compute thru pl
% in = q vector
% out = pl vector

out = pl2b(q2pl(in));
