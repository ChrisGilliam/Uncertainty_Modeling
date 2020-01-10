function [out] = pl2bel(in)

out = pl2b(in);
if out(1) < 1
	out = out / (1-out(1));
end
out(1) = 0;
