function [out] = b2q(in)
% computing FMT from b to q. Compute thru pl
% in = b vector
% out = q vector

out = pl2q(b2pl(in));
