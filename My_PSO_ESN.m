%% ����Ⱥ�Ż�Ȩ�صĻ���״̬������
% by FuMing 2015.6.3
clear all;
clc;

%% ��������
a=0.75;
N=1000;                  %�������ڲ���Ԫ����
L=1;                    %�����Ԫ����
K=10;                  %Ƕ��ά������Ϊ������Ԫ����
xs=0.03;                %����ϡ���
trainlength=1000;
testlength=100;
M0=200;                 %��M0����ʼ�ռ�����������
M=trainlength-K;
load ESNfile ESNData;   %��������
data=ESNData;
%��������ݹ�һ��,ԭ����Ԫ��x����y=(x-min)*(0.9-0.1)/(max-min)+0.1,��ʹֵ��[0.1,0.9]֮�䣬����ٰ�y����
%data_guiyi=(data-min(data)).*(0.9-0.1)./(max(data)-min(data)+0.1);
traindata=data(1:trainlength+testlength);
testdata=data(trainlength+1:trainlength+testlength);

%% ����ѵ��������Ԥ������
matrix=[];          %ÿһ��Ϊһ��ѵ��������200������һ���������900��
for i=1:M+testlength
    temp=[];
    for j=0:K
        temp=[temp,traindata(i+j)];
    end
    matrix=[matrix;temp];
end
U=matrix(1:M,1:K);  %ÿһ��Ϊһ�����룬k�����룬��800��
y=matrix(1:M,K+1);  %ÿһ��Ϊһ�������1���������800��
Ustar=matrix(M+1:M+testlength,1:K);     %ÿһ��Ϊһ�����룬k�����룬��100��
ystar=matrix(M+1:M+testlength,K+1);     %ÿһ��Ϊһ�������1���������100��

%% ����ϡ���Ϊxs��ϡ�����W0
W0=zeros(N,N);              %�ڲ�����Ȩֵ����N*N
v1=randperm(N*N);
v2=rand(N*N*xs,1);
for k=1:N*N*xs
    W0(v1(k))=v2(k,1);      %ΪW0�����ֵ���ܸ�ֵ����ΪN*N*xs��������λ��Ϊ0
end

%% ����������״̬����X
W=W0*(1/vrho(W0))*a;        %��W0*�������ӱ任�õ�W;vrho()���װ뾶����
Wfb=rand(N,L);              %��������ڲ�״̬������Ȩ�أ���������״̬��ķ���ֵ
Win=rand(N,K);              %��������ڲ�״̬������Ȩ��

X(:,1)=tanh(Win*U(1,:)');   %������״̬����������*����Ȩ�صõ�;ʹ��tanh()˫�����к�������
yn=wgn(N,1,0.0001);     %������˹������; wgn(m,n,p) ����һ��m��n�еĸ�˹�������ľ���p��dBWΪ��λָ�����������ǿ��
%�������´�����״̬����X
for i=2:M
    X(:,i)=tanh(Win*U(i,:)'+W*X(:,i-1)+Wfb*y(i-1,:)'+yn);
end

%% ESN��ѵ�����̾��Ǹ��ݸ�����ѵ������ȷ��ϵ���������Ȩ����Wout�Ĺ���
%��M0���ռ�������������������������
S=[];
D=[];
for i=M0:M
    Xnew=[X(:,i);U(i,:)']';
    S=[S;Xnew];
    D=[D;y(i,:)];
end
Wout=(pinv(S'*S)*S'*D)';

Y=(S'*Wout);
%% ʹ��PSO���Wout
G =200;   %��������
n = N+K;   %����ά��
m = 20;   %��Ⱥ��ģ
w = 0.1;  %�㷨����
c1 = 2;   %�㷨����
c2 = 2;   %�㷨����

%��ʼ����Ⱥpop
pop(m,n)=0;
for i=1:m
    pop(i,:) = Wout+0.01*rand(1,n);
end
%��ʼ�������ٶ�
V = 0.1*rands(m,n);

%���ݳ�ʼ������Ⱥ�������û�,�ҳ�Ⱥ�����ź͸�������
for s = 1:m
    indivi = pop(s,:);    %�������
    for i=M0:M
        y_f(i,:)=indivi*[X(:,i);U(i,:)'];
    end
    Er=y_f(M0:M,:)-y(M0:M,:);
    Er;
    E=(sum(Er.^2));
    E;
    Error(s) = E;
end

[OderEr,IndexEr] = sort(Error);    %������������
Error;
Errorleast = OderEr(1);    %�����С���
for i = 1:m
    if Errorleast == Error(i)
        gbest = pop(i,:);   %�ҳ���С����Ӧ�ĸ��弫ֵgbest(ָ���Ǹ���)
        break;
    end
end
ibest = pop;   %�ѳ�ʼ������Ⱥ��ΪȺ�弫ֵibest

%ѭ����ʼ
for kg = 1:G
    kg
    for s = 1:m;
        %������4%�ı������
        for j = 1:n
            for i = 1:m
                if rand(1)<0.04
                    pop(i,j) = rands(1);  %�Ը���pop(i,j)���б���
                end
            end
        end
        %r1,r2Ϊ����Ⱥ�㷨����
        r1 = rand(1);
        r2 = rand(1);
        
        % �ٶȸ���
        V(s,:) = w*V(s,:) + c1*r1*(ibest(s,:)-pop(s,:)) + c2*r2*(gbest-pop(s,:));
        %�������
        pop(s,:) = pop(s,:) + 0.3*V(s,:);
        
        %����º��ÿ���������ɿ�����Ӧ��ֵ
        for i=M0:M
            y_ff(i,:)=pop(s,:)*[X(:,i);U(i,:)'];
        end
        Er=y_ff(M0:M,:)-y(M0:M,:);
        Er;
        E=(sum(Er.^2));
        error(s) = E;
        %������Ӧ��ֵ�Ը������ź�Ⱥ�����Ž��и���
        if error(s)<Error(s)
            ibest(s,:) = pop(s,:);
            Error(s) = error(s);
        end
        if error(s)<Errorleast
            gbest = pop(s,:);
            Errorleast = error(s);
        end
    end
    
    Best(kg) = Errorleast;
end
figure(1)
plot(Best,'black');
legend('����Ⱥ�㷨ÿһ����С���');
xlabel('��������');
ylabel('Ԥ�����');
%%
%ESNԤ��
X(:,M+1)=tanh(Win*Ustar(1,:)'+W*X(:,M)+Wfb*y(M,:)'+yn);
yuce(1,:)=Wout*[X(:,M+1);Ustar(1,:)'];
for i=2:testlength
    X(:,M+i)=tanh(Win*Ustar(i,:)'+W*X(:,M+i-1)+Wfb*y(i-1,:)'+yn);
    yuce(i,:)=Wout*[X(:,M+i);Ustar(i,:)'];
end
%ESN-PSOԤ��
X(:,M+1)=tanh(Win*Ustar(1,:)'+W*X(:,M)+Wfb*y(M,:)'+yn);
yuce_pso(1,:)=gbest*[X(:,M+1);Ustar(1,:)'];
for i=2:testlength
    X(:,M+i)=tanh(Win*Ustar(i,:)'+W*X(:,M+i-1)+Wfb*y(i-1,:)'+yn);
    yuce_pso(i,:)=gbest*[X(:,M+i);Ustar(i,:)'];
end

figure(2)
plot(ystar,'b');
hold on;
plot(yuce,'rx-.');
hold on;
plot(yuce_pso,'g');
title('�����㷨��Ӧ��ϵͳԤ�����')
legend('��ʵֵ','ESNԤ��ֵ','ESN-PSOԤ��ֵ');









