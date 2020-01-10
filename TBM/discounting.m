function [m_star] = discounting(alpha,m)
%
m_star = alpha * m;
m_star(end) = 1 - alpha + m(end);

end

