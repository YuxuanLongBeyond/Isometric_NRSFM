function [pidx,pidx1,pidx2]=InvEdgeList(edges,i)

pidx1=find(edges(:,1)==i);
pidx2=find(edges(:,2)==i);
pidx=[pidx1;pidx2];

end