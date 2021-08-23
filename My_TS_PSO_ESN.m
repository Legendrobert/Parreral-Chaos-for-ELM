%% ������������Ⱥ�Ż�Ȩ�صĻ���״̬������
% by FuMing 2015.6.5
clear all;
clc;

%% ��������
a=0.75;
N=1000;                  %�������ڲ���Ԫ����
L=1;                    %�����Ԫ����
K=10;                  %Ƕ��ά������Ϊ������Ԫ����
xs=0.03;                %����ϡ���
trainlength=1152;
testlength=288;
M0=200;                 %��M0����ʼ�ռ�����������
M=trainlength-K;
load WindLoad WindLoadTrain;   %����ѵ������
load WindLoad WindLoadForecast;   %����Ԥ������
traindata=WindLoadTrain;
testdata=WindLoadForecast;
%��������ݹ�һ��,ԭ����Ԫ��x����y=(x-min)*(0.9-0.1)/(max-min)+0.1,��ʹֵ��[0.1,0.9]֮�䣬����ٰ�y����
%data_guiyi=(data-min(data)).*(0.9-0.1)./(max(data)-min(data)+0.1);

%% ����ѵ��������Ԥ������
trainmatrix=[];          %����ѵ�����飬ÿһ��Ϊһ��ѵ��������K������1���������M��
for i=1:M
    temp=[];
    for j=0:K
        temp=[temp,traindata(i+j)];
    end
    trainmatrix=[trainmatrix;temp];
end
testmatrix=[];          %����Ԥ�����飬ÿһ��Ϊһ��ѵ��������K������1���������testlength-K��
for i=1:testlength-K
    temp=[];
    for j=0:K
        temp=[temp,testdata(i+j)];
    end
    testmatrix=[testmatrix;temp];
end
U=trainmatrix(1:M,1:K);  %ÿһ��Ϊһ�����룬k�����룬��M��
y=trainmatrix(1:M,K+1);  %ÿһ��Ϊһ�������1���������M��
Ustar=testmatrix(1:testlength-K,1:K);     %ÿһ��Ϊһ�����룬k�����룬��testlength-K��
ystar=testmatrix(1:testlength-K,K+1);     %ÿһ��Ϊһ�������1���������testlength-K��

%% ����ϡ���Ϊxs��ϡ�����W0
W0=zeros(N,N);              %�ڲ�����Ȩֵ����N*N
v1=randperm(N*N);
v2=rand(N*N*xs,1);
for k=1:N*N*xs
    W0(v1(k))=v2(k,1);      %ΪW0�����ֵ���ܸ�ֵ����ΪN*N*xs��������λ��Ϊ0
end

%% ����������״̬����X
dd=vrho(W0);            %vrho()���װ뾶����
W=W0*(1/dd)*a;          %��W0*�������ӱ任�õ�W
Wfb=rand(N,L);              %��������ڲ�״̬������Ȩ��
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

%% ʹ��PSO���Wout
G =200;     %��������
n = N+K;    %����ά��
m = 20;     %��Ⱥ��ģ
w = 0.1;    %�㷨����
c1 = 2;     %�㷨����
c2 = 2;     %�㷨����
p=0.04;     %ÿ��������

limitL=0.9; %��Ӧ�ȷ����ֵ����
limitR=1.1; %��Ӧ�ȷ����ֵ����

TL=[];      %���ɱ�
du=100;     %���ڱ���Ľ��ɱ���
NS=200;     %���ڱ����һ�������е�Ԫ�ظ���

%��ʼ����Ⱥpop
pop(m,n)=0;
for i=1:m
    pop(i,:) = Wout+0.01*rand(1,n);
end
%��ʼ�������ٶ�
V = 0.1*rands(m,n);

%���ݳ�ʼ������Ⱥ�������û�,�ҳ�Ⱥ�����ź͸�������
for s = 1:m
    indivi = pop(1,:);    %�������
    for i=M0:M
        y_f(i,:)=indivi*[X(:,i);U(i,:)'];
    end
    Er=y_f(M0:M,:)-y(M0:M,:);
    Er;
    E=(sum(Er.^2));
    E;
    Error(s) = E;
end
ErrorBefore=Error;
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
    if kg==1
        for s = 1:m;
            %������p�ı������
            for j = 1:n
                for i = 1:m
                    if rand(1)<p
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
    else
        if(var(Error)/var(ErrorBefore)>limitL&&var(Error)/var(ErrorBefore)<limitR)
            %�������½���TS�Ż�
            kw=sum(ErrorBefore)/m;   %������Ӧ��ˮƽȡ�������ڸ�����Ⱥ����Ӧ��ƽ��ֵ
            ErrorBefore=Error;       %���ϴ�Error��ֵ��ErrorBefore����Ϊ�������
            for s = 1:m;
                %������p�ı������
                for j = 1:n
                    for i = 1:m
                        if rand(1)<p
                            pop(i,j) = rands(1);  %�Ը���pop(i,j)���б���
                        end
                    end
                end
                xbest=pop(s,:);                        %��Ѹ����ʼ��
                fbest=ErrorBefore(s);                  %�����Ӧ�ȳ�ʼ��
                for tw0=1:NS            %NS����һ�������ڵ���������
                    if(mod(tw0,200)==0&&mod(s,20)==0)
                       m
                       tw0
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
                    
                    %���ƶ���ĸ������Ӧ�������ˮƽ�Ƚ�
                    if E<kw%�ƶ�����Ӧ�����ڿ���ˮƽ
                        xbest=pop(s,:);
                        pp0=xbest;
                        fbest=E;
                        kw=E;
                        sz=size(TL,2);
                        if sz>du %���������У����ɱ�ĳ����ǹ̶��ģ����������ȼ������������������Ͱ�ǰ����ȥ���������о���Ϊ����
                            dd=sz-du;
                            TL(:,1:dd)=[];
                        end
                        %�����ж�����E<kwʱ��E�Ƿ��ڽ��ɱ�Χ��,�����ڣ�������ɱ�
                        s3=[];
                        for tz=1:size(TL,2)
                            if E<=TL(1,tz)-10^(-4) | E>=TL(1,tz)+10^(-4)
                                s3(tz)=0;
                            else
                                s3(tz)=1;
                            end
                        end
                        s3;
                        if sum(s3)==0   %˵��aaaa(g1,:)���ڽ��ɱ���
                            TL=[TL E];
                        end  %���ˣ��ж����
                    else%�ƶ�����Ӧ�����ڿ���ˮƽ
                        s3=[];%s3ֻ�����ڴ�������TL��ÿһ�����ıȽϣ����Ե�pp1�ĵ�һ������TL��ÿһ�����Ƚ��꣬��һ��ѭ����󣬰����ÿգ���ʼ��һѭ��
                        for tz=1:size(TL,2)
                            if E<=TL(1,tz)-10^(-4) | E>=TL(1,tz)+10^(-4)
                                s3(tz)=0;
                            else
                                s3(tz)=1;
                            end
                        end
                        s3;
                        if sum(s3)==0
                            xbest=pop(s,:);
                            pp0=xbest;
                            fbest=E;
                            TL=[TL E];
                        else
                            pp0=xbest;
                        end
                    end
                    
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
            end
        else
            ErrorBefore=Error;       %���ϴ�Error��ֵ��ErrorBefore����Ϊ�������
            for s = 1:m;
                %������p�ı������
                for j = 1:n
                    for i = 1:m
                        if rand(1)<p
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
        end
    end
    Best(kg) = Errorleast;
end
figure(1)
plot(Best,'black');
legend('������������Ⱥ�㷨ÿһ����С���');
xlabel('��������');
ylabel('Ԥ�����');
%%
%ESN-PSOԤ��
X(:,M+1)=tanh(Win*Ustar(1,:)'+W*X(:,M)+Wfb*y(M,:)'+yn);
yuce_pso(1,:)=gbest*[X(:,M+1);Ustar(1,:)'];
for i=2:testlength-K
    X(:,M+i)=tanh(Win*Ustar(i,:)'+W*X(:,M+i-1)+Wfb*y(i-1,:)'+yn);
    yuce_pso(i,:)=gbest*[X(:,M+i);Ustar(i,:)'];
end

figure(1)
plot(ystar,'b');
hold on;
plot(yuce_pso,'rx-.');
title('ϵͳԤ�����')
legend('��ʵֵ','TS-PSO-ESNԤ��ֵ');
yuce_TS_PSO_ESN=yuce_pso';
save result_TS_PSO_ESN yuce_TS_PSO_ESN;









