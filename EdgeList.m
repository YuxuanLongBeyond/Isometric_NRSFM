function [edges]=EdgeList(n,m,dist)
edges=[];
for i=1:n
    for j=1:m
        for k=1:dist
            if(j+k<=m)
            idx1=(j-1)*n+i;
            idx2=(j-1+k)*n+i;
            edges=[edges;[idx1,idx2]];
            end
            if(i+k<=n)
            idx1=(j-1)*n+i;
            idx2=(j-1)*n+i+k;
            edges=[edges;[idx1,idx2]];
            end
            if(i+k<=n & j+k<=m)
            idx1=(j-1)*n+i;
            idx2=(j-1+k)*n+i+k;
            edges=[edges;[idx1,idx2]];
            end
        end        
    end
end

p=find(edges(:,1)>n*m);

edges(p,:)=[];
p=find(edges(:,2)>m*m);
edges(p,:)=[];
end
