%% 原始数据的输入
K=10;                  %嵌入维数，即为输入神经元个数
trainlength=1152;
testlength=288;
M=trainlength-K;
load WindLoad WindLoadTrain;   %加载训练数据
load WindLoad WindLoadForecast;   %加载预测数据
traindata=WindLoadTrain;
testdata=WindLoadForecast;
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
U=trainmatrix(1:M,1:K)';  %每一行为一组输入，k个输入，共M行
y=trainmatrix(1:M,K+1)';  %每一行为一组输出，1个输出，共M行
Ustar=testmatrix(1:testlength-K,1:K)';     %每一行为一组输入，k个输入，共testlength-K行
ystar=testmatrix(1:testlength-K,K+1)'; 

%% 建立BP模型
net =newff(U,y,5);
net.trainParam.epochs = 1000;       %最大训练轮回
net.trainParam.goal = 1e2;  %均方误差
net = train(net,U,y);              %开始训练

%% 利用原始数据对BP网络仿真
yuce_BP = sim(net,Ustar);                   %用训练好的模型进行仿真
%% 图形绘制
figure(1)
plot(ystar,'b');
hold on;
plot(yuce_BP,'rx-.');
title('系统预测输出')
legend('真实值','BP预测值');

save result_BP yuce_BP;
save result_his ystar;