function m = pign2m(betP)
% Build the LC bba from the pignistic probability

n = length(betP);
two_n = 2^n;

[p,i] = sort(betP);
m = zeros(two_n,1);
last = 2^n-1;
m(last+1) = n*p(1);
for k=2:n
   ind = i(k-1);
   last = bitset(last,ind,0);
   m(last+1) = (n-k+1)*(p(k)-p(k-1));
end

   
