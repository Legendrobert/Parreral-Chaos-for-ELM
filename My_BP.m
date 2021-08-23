%% ԭʼ���ݵ�����
K=10;                  %Ƕ��ά������Ϊ������Ԫ����
trainlength=1152;
testlength=288;
M=trainlength-K;
load WindLoad WindLoadTrain;   %����ѵ������
load WindLoad WindLoadForecast;   %����Ԥ������
traindata=WindLoadTrain;
testdata=WindLoadForecast;
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
U=trainmatrix(1:M,1:K)';  %ÿһ��Ϊһ�����룬k�����룬��M��
y=trainmatrix(1:M,K+1)';  %ÿһ��Ϊһ�������1���������M��
Ustar=testmatrix(1:testlength-K,1:K)';     %ÿһ��Ϊһ�����룬k�����룬��testlength-K��
ystar=testmatrix(1:testlength-K,K+1)'; 

%% ����BPģ��
net =newff(U,y,5);
net.trainParam.epochs = 1000;       %���ѵ���ֻ�
net.trainParam.goal = 1e2;  %�������
net = train(net,U,y);              %��ʼѵ��

%% ����ԭʼ���ݶ�BP�������
yuce_BP = sim(net,Ustar);                   %��ѵ���õ�ģ�ͽ��з���
%% ͼ�λ���
figure(1)
plot(ystar,'b');
hold on;
plot(yuce_BP,'rx-.');
title('ϵͳԤ�����')
legend('��ʵֵ','BPԤ��ֵ');

save result_BP yuce_BP;
save result_his ystar;