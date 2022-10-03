%read file
function [A, b, x, M_mitgcm, M_mitgcm1, M_diagonal] = init_cg2d(cg2dpcoffdfac)

ac2d_grid=rdmds('aC2d.0000000001');
as2d_grid=rdmds('aS2d.0000000001');
aw2d_grid=rdmds('aW2d.0000000001');
aw2d=zeros(21,17);aw2d(1:20,1:16)=aw2d_grid;
as2d=zeros(21,17);as2d(1:20,1:16)=as2d_grid;
ac2d=zeros(21,17);ac2d(1:20,1:16)=ac2d_grid;

b=rdmds('cg2d_b.0000000001');
b=reshape(b,[320,1]);%右端项
x=rdmds('cg2d_x.0000000001');
x=reshape(x,[320,1]);

%A
A = zeros(20*16,20*16);
for j = 1:16
    for i =1:20
        k = (j-1)*20 +i;
        A(k,k)=ac2d(i,j);
        if k-1>=1
            A(k,k-1)=aw2d(i,j);
        end
        if k-20>=1
            A(k,k-20)=as2d(i,j);
        end
        if k+1<=20*16
            A(k,k+1)=aw2d(i+1,j);
        end
        if k+20<=20*16
            A(k,k+20)=as2d(i,j+1);
        end
    end
end

pC = zeros(21,17);pW = zeros(21,17);pS = zeros(21,17);
for j = 1:16
    for i = 1:20
        ac  = ac2d(i,j);
        if(ac==0)
            pC(i,j)=1.0;
        else
            pC(i,j)=1.0/ac;
        end
        if(i==1)
            acw=0;
        else
            acw = ac2d(i-1,j);
        end
        if(ac+acw==0)
            pW(i,j)=0.0;
        else
            pW(i,j)=-aw2d(i,j)/((cg2dpcoffdfac*(acw+ac))^2);
        end
        if(j==1)
            acs=0;
        else
            acs = ac2d(i,j-1);
        end
        if(ac+acs==0)
            pS(i,j)=0.0;
        else
            pS(i,j)=-as2d(i,j)/((cg2dpcoffdfac*(acs+ac))^2);
        end
    end
end

%使用对角线作为预处理器
M_diagonal = zeros(20*16,20*16);%preconditioner
for j = 1:16
    for i =1:20
        k = (j-1)*20 +i;
        M_diagonal(k,k)=1.0/ac2d(i,j);
    end
end

%使用mitgcm里的预处理器
M_mitgcm = zeros(20*16,20*16);
for j = 1:16
    for i =1:20
        k = (j-1)*20 +i;
        M_mitgcm(k,k)=pC(i,j);
        if k-1>=1
            M_mitgcm(k,k-1)=pW(i,j);
        end
        if k-20>=1
            M_mitgcm(k,k-20)=pS(i,j);
        end
        if k+1<=20*16
            M_mitgcm(k,k+1)=pW(i+1,j);
        end
        if k+20<=20*16
            M_mitgcm(k,k+20)=pS(i,j+1);
        end
    end
end
pC1 = zeros(21,17);pW1 = zeros(21,17);pS1 = zeros(21,17);
for j = 1:16
    for i = 1:20
        ac  = ac2d(i,j);
        if(ac==0)
            pC1(i,j)=1.0;
        else
            pC1(i,j)=1.0/ac;
        end
        if(i==1)
            acw=0;
        else
            acw = ac2d(i-1,j);
        end
        if(ac==0 || acw==0)
            pW1(i,j)=0.0;
        else
            pW1(i,j)=-aw2d(i,j)/(acw*ac);
        end
        if(j==1)
            acs=0;
        else
            acs = ac2d(i,j-1);
        end
        if(ac==0 || acs==0)
            pS1(i,j)=0.0;
        else
            pS1(i,j)=-as2d(i,j)/(acs*ac);
        end
    end
end
%使用mitgcm里的预处理器
M_mitgcm1 = zeros(20*16,20*16);
for j = 1:16
    for i =1:20
        k = (j-1)*20 +i;
        M_mitgcm1(k,k)=pC1(i,j);
        if k-1>=1
            M_mitgcm1(k,k-1)=pW1(i,j);
        end
        if k-20>=1
            M_mitgcm1(k,k-20)=pS1(i,j);
        end
        if k+1<=20*16
            M_mitgcm1(k,k+1)=pW1(i+1,j);
        end
        if k+20<=20*16
            M_mitgcm1(k,k+20)=pS1(i,j+1);
        end
    end
end
end