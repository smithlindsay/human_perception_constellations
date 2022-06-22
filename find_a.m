function [SAR] = find_a(part1,part2)

if size(part1,1)==1,
    part1=part1';
end
if size(part2,1)==1,
    part2=part2';
end
if length(part1)~=length(part2),
    disp('ERROR: partitions not of equal length')
    return
end

%Generate contingency table and calculate row/column marginals
nij=sparse(part1+1,part2+1,1);
ni=sum(nij,2);
nj=sum(nij,1);

%Identify total number of elements, n, numbers of pairs, M, and numbers of
%classified-same pairs in each partition, M1 and M2.
n=length(part1);
M=n*(n-1)/2;
M1=sum(ni.^2-ni)/2;
M2=sum(nj.^2-nj)/2;

%Pair counting types:
a=full(sum(sum(nij.^2-nij)))/2; %same in both

%Rand and Adjusted Rand indices:
meana=M1*M2/M;
SAR=(a-meana)/((M1+M2)/2-meana);