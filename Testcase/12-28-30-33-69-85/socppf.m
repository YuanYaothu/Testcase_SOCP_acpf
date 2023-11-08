function [v_true,R_t,X_t,node_del,err_relax1]=socppf(numbus,numbranch,branch,H,A,Pit,Qit)
con_b=branch(:,2:3);
con_b=con_b-1;
xhat = sdpvar(numbus-1+numbranch*2,1);
h= sum(xhat(numbranch+1:numbranch*2));

F=[ (xhat(numbranch+1)^2+xhat(2*numbranch+1)^2) <= 1*xhat(1) ];

for t=2:numbranch
    F= [ F, (xhat(numbranch+t)^2+xhat(2*numbranch+t)^2)<=(xhat(con_b(t,1))*xhat(con_b(t,2)))];
end

F=[F,[Pit(2:end,:);Qit(2:end,:)]-H*xhat(1:end)==0,xhat(1:2*numbranch)>=0];

options = sdpsettings('solver','mosek');
o=optimize( F , -h , options);
result=value(xhat);

v_true=sqrt(result(1:numbranch,1));
R_t=result(numbranch+1:2*numbranch,1);
X_t=result(2*numbranch+1:end,1);

del_1=zeros(numbranch,1);
del_2=zeros(numbranch,1);
del_true=zeros(numbranch,1);

err_relax1=zeros(numbranch,1);

for t=1:numbranch
    I = con_b(t, 1) ;
    J = con_b(t, 2) ;

    if I==0


        del_true(t,1)= atan(X_t(t,1)/R_t(t,1));

        err_relax1(t)=(R_t(t)^2+X_t(t)^2)-1*v_true(J,1)^2;

    else


        err_relax1(t)=(R_t(t)^2+X_t(t)^2)-(v_true(I,1)^2)*(v_true(J,1)^2);
        del_true(t,1)= atan(X_t(t,1)/R_t(t,1));


    end
end

A(:,1)=[];
node_del=inv(A)*del_true;

end