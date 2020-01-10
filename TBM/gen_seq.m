function [seq] = gen_seq(n,C)
% Generate all sequences of length n with C numbers
%  [seq] = gen_seq(n,C)

% Jan 2004

nseq = C^n;
%fprintf('Number of %d sequences\n',nseq);
seq = zeros(nseq,n);

for j=1:n
   w = C^(n-j);
   for k=1:nseq/w
      seq((k-1)*w+1:k*w,j) = rem((k-1),C)*ones(w,1);
   end
end

      
   