function [p] = pignistic(m)
% Pignistic transformmation
% TBM
% B. Ristic, RMIT 2018

two_n = length(m);
n = log2(two_n);


p = zeros(two_n,1);
epsilon = 1e-10;
ind = find(m > epsilon);
for i=1:length(ind)
    a = ind(i)-1;
    if a > 0
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
            p(hh(h)+1) = p(hh(h)+1) + m(a+1)/card;
        end
    end
end
if m(1) == 1,  m(1) = 1-epsilon; end
p = p / (1-m(1));
%if sum(p)< 0.9999; keyboard; end



