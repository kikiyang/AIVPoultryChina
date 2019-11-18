% complete the nodelist according to dictionary matrix

function [wnl,freqm,Adjlist] = orignodelist2adjm(onlm,dictlist)
N = length(onlm);
onnode = size(onlm,2)-1;
wnl = cell(N,1);
freqm = zeros(32);
Adjlist = cell(N,1);
for i=1:1:N
    wnlm(1,1) = onlm(i,1);
    adjm = zeros(32);
    for j=2:1:onnode
        A = dictlist{onlm(i,j),onlm(i,j+1)};
        if j==2
            wnlm(1,[(end+1):(end+length(A))]) = A;
        else
            wnlm(1,[end:(end+length(A)-1)]) = A;
        end
    end
    for j =2:1:length(wnlm)-1
        tmpm = zeros(32);
        tmpm(wnlm(1,j),wnlm(1,j+1)) = 1;
        adjm = adjm + tmpm;
    end
    wnl{i} = wnlm;
    % generate the frequency matrix and adjacency matrix
    freqm = freqm + adjm;
    Adjlist{i} = adjm;
    clear wnlm;
end