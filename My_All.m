clear all
load result_TS_PSO_ESN yuce_TS_PSO_ESN;    %����TS_PSO_ESNԤ������
load result_ESN yuce_ESN;                  %����ESNԤ������
load result_RBF yuce_RBF;                  %����RBFԤ������
load result_BP yuce_BP;                    %����BPԤ������
load result_his ystar;                     %������ʷ����

error_BP=abs(yuce_BP-ystar)./ystar;
error_RBF=abs(yuce_RBF-ystar)./ystar;
error_ESN=abs(yuce_ESN-ystar)./ystar;
error_TS_PSO_ESN=abs(yuce_TS_PSO_ESN-ystar)./ystar;

error_avg_BP=mean(error_BP)
error_avg_RBF=mean(error_RBF)
error_avg_ESN=mean(error_ESN)
error_avg_TS_PSO_ESN=mean(error_TS_PSO_ESN)

figure(1)
plot(ystar,'b');
hold on;
plot(yuce_TS_PSO_ESN,'rx-.');
title('TS_PSO_ESNԤ�����')
legend('��ʵֵ','TS-PSO-ESN');
figure(2)
plot(ystar,'b');
hold on;
plot(yuce_ESN,'rx-.');
%hold on;
%plot(yuce_RBF,'yx-.');
hold on;
plot(yuce_BP,'g');
title('�Աȷ���Ԥ�����')
legend('��ʵֵ','ESN','BP');
figure(3)
plot(error_BP,'b');
hold on;
%plot(error_RBF,'y');
%hold on;
plot(error_ESN,'g');
hold on;
plot(error_TS_PSO_ESN,'r');
title('������������')
%legend('BP','RBF','ESN','TS-PSO-ESN');
legend('BP','ESN','TS-PSO-ESN');

plot(error_BP,'b');
hold on;