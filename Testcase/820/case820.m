% SOCP for AC power flow calculation
% Test case: 820 distribution power system

clc;
clear;

bus = load('bus820.mat');
bus = bus.bus;

branch = load('branch819.mat');
branch = branch.branch;


numbus=size(bus,1);        
numbranch=size(branch,1);  

branch_ts=[find(branch(:,9)~=0),branch(branch(:,9)~=0,:)];
m_ts=size(branch_ts,1);  

branch_pt=[find(branch(:,9)==0),branch(branch(:,9)==0,:)];
branch_pt(:,10)=1;

slack_bus=bus(bus(:,2)==3,1);
PQ_bus=bus(bus(:,2)==1,1);


Y=zeros(numbus,numbus);
for i=1:numbranch
    I=branch(i,1);  
    J=branch(i,2);  
    
    r=branch(i,3); % resistance
    X=branch(i,4); % reactance
    b=branch(i,5); % total line charging susceptance

    y=1/(r+1i*X);

    branch(branch(:,9)==0,9)=1;
   
    TAP=branch(i,9);   

    Y(I,I)=Y(I,I)+(y+1i*(b/2))/TAP^2;
    Y(J,J)=Y(J,J)+(y+1i*(b/2));
    Y(I,J)=-y/TAP;
    Y(J,I)=Y(I,J);
end

for j=1:numbus
    y0=bus(j,5)+1i*bus(j,6);
    Y(j,j)=Y(j,j)+y0;
end
G=real(Y);
B=imag(Y);


Sn=100; % MW


Pit=bus(2:end,3)/Sn;
Pit(PQ_bus-1,:)=-Pit(PQ_bus-1,:);

Qit=bus(2:end,4)/Sn;
Qit=-Qit;
 

H=zeros(2*numbus,2*numbranch+numbus);

for i=1:numbus   
   H(i,i)= G(i,i);         
   H(numbus+i,i)=-B(i,i);  
end

A=zeros(numbranch,numbus);
for j=1:numbranch
    I = branch(j, 1) ;      
    J = branch(j, 2) ;      
    
    A(j,I)=1;  
    A(j,J)=-1;
    
    H(I,numbus+j)=G(I,J);             
    H(I,numbus+numbranch+j)=B(I,J);   
    H(J,numbus+j)=G(J,I);             
    H(J,numbus+numbranch+j)=-B(J,I);  
    
    H(numbus+I,numbus+j)=-B(I,J);                                                                           
    H(numbus+I,numbus+numbranch+j)= G(I,J);  
    H(numbus+J,numbus+j)=-B(J,I);            
    H(numbus+J,numbus+numbranch+j)=-G(J,I);  
    
end
     
con_b=branch(:,1:2);  
con_b=con_b-1;
xhat = sdpvar(numbus-1+numbranch*2,1);  
h= sum(xhat(numbranch+1:numbranch*2));  

F=[ (xhat(numbranch+1)^2+xhat(2*numbranch+1)^2) <= 1*xhat(1) ];

for i=2:numbranch
    if con_b(i,1)==0
        F= [ F, (xhat(numbranch+i)^2+xhat(2*numbranch+i)^2)<=(1*xhat(con_b(i,2)))];
    else
        F= [ F, (xhat(numbranch+i)^2+xhat(2*numbranch+i)^2)<=(xhat(con_b(i,1))*xhat(con_b(i,2)))];
    end
end

H_1=H;

H(1,:)=[];
H(numbus,:)=[];
H(:,1)=[];

F=[F,[Pit;Qit]-H*xhat(1:end)==0,xhat(1:2*numbranch)>=0];

% AC power flow calculation based on MOSEK

options = sdpsettings('solver','mosek','CACHESOLVERS','1');
o=optimize( F , -h , options);
result=value(xhat); 

v_true=sqrt(result(1:numbranch,1));  
R_t=result(numbranch+1:2*numbranch,1);
X_t=result(2*numbranch+1:end,1);
del_1=zeros(numbranch,1);


err_relax1=zeros(numbranch,1);

v_true=[1;v_true];

for i=1:numbranch
   I = branch(i, 1) ;
   J = branch(i, 2) ;
   
   err_relax1(i)=(R_t(i)^2+X_t(i)^2)-(v_true(I,1)^2)*(v_true(J,1)^2);
   
end

ave_err=sum(err_relax1)/size(err_relax1,1);






 