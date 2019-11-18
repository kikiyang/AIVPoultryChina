% node list to adjacency matrix
function [freqm,adjm] = nodelist2adjm(nodelistm,nnode,nbrknode)
N = size(nodelistm,1);
adjm = cell(N,1);
freqm = zeros(nnode);
for i=1:1:N
    A = zeros(nnode);
    switch nbrknode
        case 0
            for j=2:1:nnode
                A(nodelistm(i,j),nodelistm(i,j+1))=1;
            end
        case 1
            brknode = nodelistm(i,end);
            for j=2:1:nnode
                if nodelistm(i,j)~=brknode,
                    A(nodelistm(i,j),nodelistm(i,j+1))=1;
                end
            end
        case 2
            brknode(i,[1,2])=[nodelistm(i,end-1),nodelistm(i,end)];
            for j=2:1:nnode
                if nodelistm(i,j)~=brknode(i,1) && nodelistm(i,j)~=brknode(i,2),
                    A(nodelistm(i,j),nodelistm(i,j+1))=1;
                end
            end
        otherwise
            break;
    end
    adjm{i,1} = A;
    freqm = freqm + A;
end
        