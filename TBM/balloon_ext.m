function m1 = balloon_ext(m,nn1,nn2)
%
% m1 = balloon_ext(m,nn1,nn2)
%
% m - bba on frame H' subset of H
% nn1+nn2 = card(H) - card(H');
% m1 - bba on frame H
% ex: m = [0 0.4 0.2 0.4]'; m1 = balloon_ext(m,2,1);

two_n = length(m);
n = log2(two_n);  % cardinality of H'

N = n + nn1 + nn2;   		% cardinality of H
two_N = 2^N;

m1 = zeros(two_N,1);
for i=1:two_n
   ind=i-1;  % binary index in H
   ind1 = convert(ind,n,nn1,nn2);
   m1(ind1+1) = m(ind+1); 
end

%-------------------------------------------
function a = convert(b,n,nn1,nn2);
%
a = 0;
if nn1 > 0
   for i=1:nn1
      a = bitset(a,i);
   end
end
j=0;
for i=1+nn1:n+nn1
   j=j+1;
   c = bitget(b,j);
   a = bitset(a,i,c);
end
if nn2 > 0
   for i=1+nn1+n:nn1+nn2+n
      a = bitset(a,i);
   end
end