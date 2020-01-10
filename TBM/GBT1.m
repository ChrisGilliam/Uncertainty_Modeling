function [mH] = GBT1(mZ_table,measZ)
%
% Implementation of the GBT according to the paper by Delmotte & Smets
% Input: mZ_table is a bba table conditioned on Hypotheses (X) in the 
%				Mesaurement space (Z)
%        measZ = index on Z representing the measurement (e.g. 3 means z=2,
%				4 is {1,2}, etc
% Output: m (bba) conditioned on measurement (Z) in the Hypotheses space (X)
% Example:
% mm = zeros(8,3);mm(2,:)=[1 0.2 0.2];mm(4,:) =[0 0.8 0.1];mm(8,3)=0.7;
% [mH] = GBT1(mm,2);  % this is conditioned on measurement z1
% displ_bel(mH)

%addpath '.\FMT\'

epsilon = 1.0000e-010;

[two_Z_card,H_card] = size(mZ_table);

for i=1:H_card
  m = mZ_table(:,i);
  pl = m2pl(m);
  plZ_table(:,i) = pl;
end

lik = plZ_table(measZ,:);
ind = find(lik == 1);
lik(ind) = 1 - epsilon;

mH = prod(1 -lik);
for i=1:H_card
   mH = [mH mH*lik(i)/(1 - lik(i))];
end
mH = mH';


