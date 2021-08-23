%% 禁忌搜索粒子群优化权重的回声状态神经网络
% by FuMing 2015.6.5
clear all;
clc;

%% 参数设置
a=0.75;
N=1000;                  %储备池内部神经元个数
L=1;                    %输出神经元个数
K=10;                  %嵌入维数，即为输入神经元个数
xs=0.03;                %矩阵稀疏度
trainlength=1152;
testlength=288;
M0=200;                 %从M0步开始收集储备池向量
M=trainlength-K;
load WindLoad WindLoadTrain;   %加载训练数据
load WindLoad WindLoadForecast;   %加载预测数据
traindata=WindLoadTrain;
testdata=WindLoadForecast;
%下面把数据归一化,原来的元素x换成y=(x-min)*(0.9-0.1)/(max-min)+0.1,以使值在[0.1,0.9]之间，最后再把y换回
%data_guiyi=(data-min(data)).*(0.9-0.1)./(max(data)-min(data)+0.1);

%% 构造训练样本、预测样本
trainmatrix=[];          %构造训练数组，每一行为一个训练样本，K个输入1个输出，共M行
for i=1:M
    temp=[];
    for j=0:K
        temp=[temp,traindata(i+j)];
    end
    trainmatrix=[trainmatrix;temp];
end
testmatrix=[];          %构造预测数组，每一行为一个训练样本，K个输入1个输出，共testlength-K行
for i=1:testlength-K
    temp=[];
    for j=0:K
        temp=[temp,testdata(i+j)];
    end
    testmatrix=[testmatrix;temp];
end
U=trainmatrix(1:M,1:K);  %每一行为一组输入，k个输入，共M行
y=trainmatrix(1:M,K+1);  %每一行为一组输出，1个输出，共M行
Ustar=testmatrix(1:testlength-K,1:K);     %每一行为一组输入，k个输入，共testlength-K行
ystar=testmatrix(1:testlength-K,K+1);     %每一行为一组输出，1个输出，共testlength-K行

%% 定义稀疏度为xs的稀疏矩阵W0
W0=zeros(N,N);              %内部连接权值矩阵N*N
v1=randperm(N*N);
v2=rand(N*N*xs,1);
for k=1:N*N*xs
    W0(v1(k))=v2(k,1);      %为W0随机赋值，总赋值个数为N*N*xs个，其余位置为0
end

%% 构建储备池状态变量X
dd=vrho(W0);            %vrho()求谱半径函数
W=W0*(1/dd)*a;          %将W0*缩放因子变换得到W
Wfb=rand(N,L);              %输出层与内部状态量连接权重
Win=rand(N,K);              %输入层与内部状态量连接权重

X(:,1)=tanh(Win*U(1,:)');   %储备池状态变量；输入*连接权重得到;使用tanh()双曲正切函数处理
yn=wgn(N,1,0.0001);     %产生高斯白噪声; wgn(m,n,p) 产生一个m行n列的高斯白噪声的矩阵，p以dBW为单位指定输出噪声的强度
%样本更新储备池状态向量X
for i=2:M
    X(:,i)=tanh(Win*U(i,:)'+W*X(:,i-1)+Wfb*y(i-1,:)'+yn);
end

%% ESN的训练过程就是根据给定的训练样本确定系数输出连接权矩阵Wout的过程
%从M0步收集储备池向量和样本输入向量
S=[];
D=[];
for i=M0:M
    Xnew=[X(:,i);U(i,:)']';
    S=[S;Xnew];
    D=[D;y(i,:)];
end
Wout=(pinv(S'*S)*S'*D)';

%% 使用PSO获得Wout
G =200;     %迭代次数
n = N+K;    %粒子维数
m = 20;     %种群规模
w = 0.1;    %算法参数
c1 = 2;     %算法参数
c2 = 2;     %算法参数
p=0.04;     %每代变异率

limitL=0.9; %适应度方差比值下限
limitR=1.1; %适应度方差比值上限

TL=[];      %禁忌表
du=100;     %用于变异的禁忌表长度
NS=200;     %用于变异的一个邻域中的元素个数

%初始化种群pop
pop(m,n)=0;
for i=1:m
    pop(i,:) = Wout+0.01*rand(1,n);
end
%初始化粒子速度
V = 0.1*rands(m,n);

%根据初始化的种群计算个体好坏,找出群体最优和个体最优
for s = 1:m
    indivi = pop(1,:);    %抽出个体
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
[OderEr,IndexEr] = sort(Error);    %对误差进行排序
Error;
Errorleast = OderEr(1);    %求出最小误差
for i = 1:m
    if Errorleast == Error(i)
        gbest = pop(i,:);   %找出最小误差对应的个体极值gbest(指的是个体)
        break;
    end
end
ibest = pop;   %把初始化的种群作为群体极值ibest

%循环开始
for kg = 1:G
    kg
    if kg==1
        for s = 1:m;
            %个体有p的变异概率
            for j = 1:n
                for i = 1:m
                    if rand(1)<p
                        pop(i,j) = rands(1);  %对个体pop(i,j)进行变异
                    end
                end
            end
            %r1,r2为粒子群算法参数
            r1 = rand(1);
            r2 = rand(1);
            
            % 速度更新
            V(s,:) = w*V(s,:) + c1*r1*(ibest(s,:)-pop(s,:)) + c2*r2*(gbest-pop(s,:));
            %个体更新
            pop(s,:) = pop(s,:) + 0.3*V(s,:);
            
            %求更新后的每个个体误差，可看成适应度值
            for i=M0:M
                y_ff(i,:)=pop(s,:)*[X(:,i);U(i,:)'];
            end
            Er=y_ff(M0:M,:)-y(M0:M,:);
            Er;
            E=(sum(Er.^2));
            error(s) = E;
            
            %根据适应度值对个体最优和群体最优进行更新
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
            %此条件下进行TS优化
            kw=sum(ErrorBefore)/m;   %渴望适应度水平取个体所在父代的群体适应度平均值
            ErrorBefore=Error;       %把上次Error赋值给ErrorBefore，作为父代误差
            for s = 1:m;
                %个体有p的变异概率
                for j = 1:n
                    for i = 1:m
                        if rand(1)<p
                            pop(i,j) = rands(1);  %对个体pop(i,j)进行变异
                        end
                    end
                end
                xbest=pop(s,:);                        %最佳个体初始化
                fbest=ErrorBefore(s);                  %最佳适应度初始化
                for tw0=1:NS            %NS是在一个邻域内的搜索次数
                    if(mod(tw0,200)==0&&mod(s,20)==0)
                       m
                       tw0
                    end
                    %r1,r2为粒子群算法参数
                    r1 = rand(1);
                    r2 = rand(1);
                    
                    % 速度更新
                    V(s,:) = w*V(s,:) + c1*r1*(ibest(s,:)-pop(s,:)) + c2*r2*(gbest-pop(s,:));
                    %个体更新
                    pop(s,:) = pop(s,:) + 0.3*V(s,:);
                    
                    %求更新后的每个个体误差，可看成适应度值
                    for i=M0:M
                        y_ff(i,:)=pop(s,:)*[X(:,i);U(i,:)'];
                    end
                    Er=y_ff(M0:M,:)-y(M0:M,:);
                    Er;
                    E=(sum(Er.^2));
                    error(s) = E;
                    
                    %把移动后的个体的适应度与渴望水平比较
                    if E<kw%移动后适应度优于渴望水平
                        xbest=pop(s,:);
                        pp0=xbest;
                        fbest=E;
                        kw=E;
                        sz=size(TL,2);
                        if sz>du %整个迭代中，禁忌表的长度是固定的，所以若长度即列数超过，超过几就把前几列去掉，操作中就设为空阵
                            dd=sz-du;
                            TL(:,1:dd)=[];
                        end
                        %下面判断满足E<kw时，E是否在禁忌表范围内,若不在，存入禁忌表
                        s3=[];
                        for tz=1:size(TL,2)
                            if E<=TL(1,tz)-10^(-4) | E>=TL(1,tz)+10^(-4)
                                s3(tz)=0;
                            else
                                s3(tz)=1;
                            end
                        end
                        s3;
                        if sum(s3)==0   %说明aaaa(g1,:)不在禁忌表内
                            TL=[TL E];
                        end  %到此，判断完毕
                    else%移动后适应度劣于渴望水平
                        s3=[];%s3只是用于存放相对于TL的每一个数的比较，所以当pp1的第一个数与TL的每一个数比较完，即一次循环完后，把它置空，开始下一循环
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
                    
                    %根据适应度值对个体最优和群体最优进行更新
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
            ErrorBefore=Error;       %把上次Error赋值给ErrorBefore，作为父代误差
            for s = 1:m;
                %个体有p的变异概率
                for j = 1:n
                    for i = 1:m
                        if rand(1)<p
                            pop(i,j) = rands(1);  %对个体pop(i,j)进行变异
                        end
                    end
                end
                %r1,r2为粒子群算法参数
                r1 = rand(1);
                r2 = rand(1);
                
                % 速度更新
                V(s,:) = w*V(s,:) + c1*r1*(ibest(s,:)-pop(s,:)) + c2*r2*(gbest-pop(s,:));
                %个体更新
                pop(s,:) = pop(s,:) + 0.3*V(s,:);
                
                %求更新后的每个个体误差，可看成适应度值
                for i=M0:M
                    y_ff(i,:)=pop(s,:)*[X(:,i);U(i,:)'];
                end
                Er=y_ff(M0:M,:)-y(M0:M,:);
                Er;
                E=(sum(Er.^2));
                error(s) = E;
                
                %根据适应度值对个体最优和群体最优进行更新
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
legend('紧急搜索粒子群算法每一代最小误差');
xlabel('进化次数');
ylabel('预测误差');
%%
%ESN-PSO预测
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
title('系统预测输出')
legend('真实值','TS-PSO-ESN预测值');
yuce_TS_PSO_ESN=yuce_pso';
save result_TS_PSO_ESN yuce_TS_PSO_ESN;









