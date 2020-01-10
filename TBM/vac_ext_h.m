function m_ext = vac_ext_h(m,nn);
% VACOUOUS EXTENSION on XxY 
% m_ext = vac_ext_h(m,nn);
% m - input bba on frame X={x1,x2,..x_n}; m is defined on 2^n
% you want to extend to (XxY), where Y = {y1,...,y_nn}

two_n = length(m);
n = log2(two_n);  % cardinality of X

N = n*nn;
two_N = 2^N;

m_ext = zeros(two_N,1);
for i=1:two_n
   if m(i) > 0
        ind=i-1;  % binary index in X
        ind_ext = convert(ind,n,nn);
        m_ext(ind_ext+1) = m(ind+1);
    end
end


%-------------------------------------------
function a = convert(b,n,nn);
%
a=0;
for i=1:nn
    c = bitshift(b,(i-1)*n);
    a = bitor(a,c);
end
