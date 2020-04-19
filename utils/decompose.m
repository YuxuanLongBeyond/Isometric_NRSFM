function [R, t, n, d] = decompose(H,x,y) 

[U D V]= svd(H); 
 
d1=D(1,1); 
d2=D(2,2); 
d3=D(3,3); 
s=det(U)*det(V); 
 
R1 = eye(3) ; t1 = zeros(3,1) ; n1 = zeros(3, 1) ; 
R2 = eye(3) ; t2 = zeros(3,1) ; n2 = zeros(3, 1) ; 
R3 = eye(3) ; t3 = zeros(3,1) ; n3 = zeros(3, 1) ; 
R4 = eye(3) ; t4 = zeros(3,1) ; n4 = zeros(3, 1) ; 
 
tmp = (H(3,1)*x + H(3,2)*y + H(3,3))/(d2*s); 
if(tmp > 0) 
    d = d2; 
else     
    d = -d2; 
end 
 
ep=0.0001; 
if (abs(d1-d2) < ep) 
    if (abs(d2-d3) < ep) 
        cas = 4; 
    else 
        cas = 2; 
    end 
else 
    if (abs(d2-d3) < ep) 
        cas = 3; 
    else 
        cas = 1; 
    end 
end 
%  cas 
%   
%   if(d<0) 
%       cas 
%       pause; 
%   end 
 
 
x1=sqrt((d1*d1 - d2*d2)/(d1*d1 - d3*d3)); 
x2=0; 
x3=sqrt((d2*d2 - d3*d3)/(d1*d1 - d3*d3)); 
 
norm1ok=0; 
norm2ok=0; 
norm3ok=0; 
norm4ok=0; 
 
if(cas==1) 
    n1=[x1 x2 x3]'; 
    N1=V*n1; 
    tmp = (N1(1)*x + N1(2)*y + N1(3))/(d*s); 
    norm1ok = (tmp>0); 
     
    n2=[x1 x2 -x3]'; 
    N2=V*n2; 
    tmp = (N2(1)*x + N2(2)*y + N2(3))/(d*s); 
    norm2ok = (tmp>0); 
     
    n3=[-x1 x2 x3]'; 
    N3=V*n3; 
    tmp = (N3(1)*x + N3(2)*y + N3(3))/(d*s); 
    norm3ok = (tmp>0); 
     
    n4=[-x1 x2 -x3]'; 
    N4=V*n4; 
    tmp = (N4(1)*x + N4(2)*y + N4(3))/(d*s); 
    norm4ok = (tmp>0); 
end 
if(cas==2)     
    %x3 = +-1 
    n1=[0 0 1]'; 
    N1 = V*n1; 
    tmp = (N1(1)*x + N1(2)*y + N1(3))/(d*s); 
    norm1ok = (tmp>0); 
     
    n2=[0 0 -1]'; 
    N2 =V*n2; 
    tmp = (N2(1)*x + N2(2)*y + N2(3))/(d*s); 
    norm2ok = (tmp>0); 
end 
if(cas==3) 
    %x1 = +-1 
    n1=[1 0 0]'; 
    N1 =V*n1; 
    tmp = (N1(1)*x + N1(2)*y + N1(3))/(d*s); 
    norm1ok = (tmp>0); 
     
    n2=[-1 0 0]'; 
    N2 =V*n2; 
    tmp = (N2(1)*x + N2(2)*y + N2(3))/(d*s); 
    norm2ok = (tmp>0); 
end 
     
if(d>0) 
    if (cas == 1) 
        sinth= sqrt((d1*d1 - d2*d2)*(d2*d2 - d3*d3))/((d1+d3)*d2); 
        costh= (d2*d2 + d1*d3)/((d1+d3)*d2); 
        % if (norm1ok)                 
            R1=[costh 0 -sinth;0 1 0; sinth 0 costh];            
            t1=(d1-d3) * [n1(1) 0 -n1(3)]';  
        %end 
        %if (norm2ok) 
            R2=[costh 0 sinth;0 1 0; -sinth 0 costh];             
            t2=(d1-d3)*[n2(1) 0 -n2(3)]'; 
        %end 
        %if (norm3ok) 
            R3=[costh 0 sinth;0 1 0; -sinth 0 costh];             
            t3=(d1-d3)*[n3(1) 0 -n3(3)]'; 
        %end 
        %if (norm4ok) 
            R4=[costh 0 -sinth;0 1 0; sinth 0 costh];             
            t4=(d1-d3)*[n4(1) 0 -n4(3)]'; 
        %end 
    end 
    if(cas == 2) 
        %if (norm1ok) 
            R1=eye(3); 
            t1=(d3-d1)*n1; 
        %end 
        %if(norm2ok) 
            R2=eye(3); 
            t2=(d3-d1)*n2; 
        %end 
    end 
    if (cas == 3) 
        %if (norm1ok) 
            R1=eye(3); 
            t1=(d1-d2)*n1; 
        %end 
        %if (norm2ok) 
            R2=eye(3); 
            t2=(d1-d2)*n2; 
        %end 
    end 
    if (cas == 4) 
        R1=eye(3); 
        t1=[0 0 0]'; 
    end     
end 
 
if(d<0) 
    if (cas == 1)         
        sinph=sqrt( (d1*d1 - d2*d2) * (d2*d2 - d3*d3) ) / ( (d1-d3) *d2); 
        cosph=(d1*d3 - d2*d2)/( (d1-d3) *d2); 
        %if (norm1ok)             
            R1=[cosph 0 sinph;0 -1 0; sinph 0 -cosph];             
            t1=(d1+d3)*n1; 
        %end 
        %if (norm2ok) 
            R2=[cosph 0 -sinph;0 -1 0; -sinph 0 -cosph];             
            t2=(d1+d3)*n2; 
        %end 
        %if (norm3ok) 
            R3=[cosph 0 -sinph;0 -1 0; -sinph 0 -cosph];             
            t3=(d1+d3)*n3; 
        %end 
        %if (norm4ok) 
            R4=[cosph 0 sinph;0 -1 0; sinph 0 -cosph];             
            t4=(d1+d3)*n4; 
        %end 
    end 
    if (cas == 2) 
        %if (norm1ok) 
            R1=eye(3); 
            R1(1,1)=-1; 
            R1(2,2)=-1;             
            t1=(d1+d3)*n1; 
        %end 
        %if (norm2ok) 
            R2=eye(3); 
            R2(1,1)=-1; 
            R2(2,2)=-1;             
            t2=(d1+d3)*n2; 
        %end 
    end 
    if(cas == 3) 
        %if(norm1ok) 
            R1=eye(3); 
            R1(2,2)=-1; 
            R1(3,3)=-1;             
            t1=(d1+d3)*n1; 
        %end 
        %if (norm2ok) 
            R2=eye(3); 
            R2(2,2)=-1; 
            R2(3,3)=-1;             
            t2=(d1+d3)*n2; 
        %end 
    end 
    % case 4 is complex to code ... :) 
end 
 
% convert solutions from R',t',n' to R,t,n 
i=1; 
if (d>0 | cas ~= 4) 
    %if (norm1ok) 
        R{i}=s*U*R1*V'; 
        t{i}=(U*t1)/(d*s); 
        n{i}=V*n1 ;       
        i=i+1; 
         
    %end 
    %if (norm2ok) 
        R{i}=s*U*R2*V'; 
        t{i}=(U*t2)/(d*s); 
        n{i}=V*n2; 
        i=i+1; 
         
    %end 
    %if (norm3ok) 
        R{i}=s*U*R3*V'; 
        t{i}=(U*t3)/(d*s); 
        n{i}=V*n3; 
        i=i+1; 
         
    %end 
    %if (norm4ok) 
        R{i}=s*U*R4*V'; 
        t{i}=(U*t4)/(d*s); 
        n{i}=V*n4; 
        i=i+1;         
    %end 
end 
 
d=d*s;
 
% verify if correct 
for i=1:length(R) 
     Hc=R{i}+t{i}*n{i}'; 
     Hc=Hc/Hc(3,3); 
     errorInHdecomposition=H-Hc; 
%     R{i}=R{i}/R{i}(3,3); 
     %t{i}=t{i}/norm(t{i}); 
 end 
 
 
 
 
 
 
 
    
 
 
 
 
 
 
 
