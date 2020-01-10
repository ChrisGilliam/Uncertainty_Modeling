function [p,bc] = voorbraak(m);
% Voorbraak transformmation

two_n = length(m);
n = log2(two_n);

% Compute bayesian constant first
bc = 0;
epsilon = 0.0000001;
ind = find(m > epsilon);
for i=1:length(ind)
    a = ind(i) - 1;
    % compute cardinality
    card = 0;
    for b=1:n
        h = 2^(b-1);
        if bitand(a,h)
           card = card + 1;
        end
    end
    bc = bc + m(a+1)*card; 
end

p = zeros(two_n,1);
for i=1:length(ind)
   a = ind(i)-1;
   card = 0;
   hh = [];
   for b=1:n
       h = 2^(b-1);
       if bitand(a,h)
           card = card + 1;
           hh(card)=h;
       end
   end
   for h=1:length(hh)
        p(hh(h)+1) = p(hh(h)+1) + m(a+1);
   end
end
p = p/bc;