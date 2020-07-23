function res = solve_normal(eq)

mpol q1
mpol q2

    
pol = [q1^4,q1^3*q2,q1^2*q2^2,q1^1*q2^3,q2^4,q1^3,q2^3,q1^2,q2^2,q1,q2,1]*eq;
P = msdp(min(pol));
[status,obj,M] = msol(P);
if status ==1
    res = double([q1 q2]);
else
    res = [0,0];
end
