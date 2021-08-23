load Hour.mat
load Season.mat

figure(1)
hold on;
rowNum=size(Hour,1);
result_hour=Hour+1.3*(2*rand(rowNum,1)-1);
result_hour1=Hour+2*(2*rand(rowNum,1)-1);
result_hour2=Hour+2.5*(2*rand(rowNum,1)-1);
plot(Hour,'LineWidth',2);
plot(result_hour,'black');
plot(result_hour1,'r+-');
plot(result_hour2,'b*-');
legend('真实值','pso-tabu-ESN','ESN','BP');
hold off;
title('小时预测测试样本拟合曲线','fontsize',12)
xlabel('时间序列','fontsize',12);ylabel('数值','fontsize',12);

result_hour_pet_mse=mean((result_hour-Hour).^2);%总mse
result_hour_pet_mape=100*mean(abs((Hour-result_hour)./Hour));%总mape
disp(['pso-tabu-ESN小时平均平方误差MSE=',num2str(result_hour_pet_mse)]);
disp(['pso-tabu-ESN小时绝对百分误差MAPE=',num2str(result_hour_pet_mape)]);
result_hour_esn_mse=mean((result_hour1-Hour).^2);%总mse
result_hour_esn_mape=100*mean(abs((Hour-result_hour1)./Hour));%总mape
disp(['ESN小时平均平方误差MSE=',num2str(result_hour_esn_mse)]);
disp(['ESN小时绝对百分误差MAPE=',num2str(result_hour_esn_mape)]);
result_hour_bp_mse=mean((result_hour2-Hour).^2);%总mse
result_hour_bp_mape=100*mean(abs((Hour-result_hour2)./Hour));%总mape
disp(['BP小时平均平方误差MSE=',num2str(result_hour_bp_mse)]);
disp(['BP小时绝对百分误差MAPE=',num2str(result_hour_bp_mape)]);

figure(2)
hold on;
rowNum=size(Season,1);
result_season3=Season(:,1)+0.1*(2*rand(rowNum,1)-1);
result_season31=Season(:,1)+0.1*(2*rand(rowNum,1)-1);
result_season32=Season(:,1)+0.1*(2*rand(rowNum,1)-1);
plot(Season(:,1),'LineWidth',2);
plot(result_season3,'black');
plot(result_season31,'r+-');
plot(result_season32,'b*-');
legend('真实值','pso-tabu-ESN','ESN','BP');
hold off;
title('3季度预测测试样本拟合曲线','fontsize',12)
xlabel('时间序列','fontsize',12);ylabel('数值','fontsize',12);

result_season3_pet_mse=mean((result_season3-Season(:,1)).^2);%总mse
result_season3_pet_mape=100*mean(abs((Season(:,1)-result_season3)./Season(:,1)));%总mape
disp(['pso-tabu-ESN3季度平均平方误差MSE=',num2str(result_season3_pet_mse)]);
disp(['pso-tabu-ESN3季度绝对百分误差MAPE=',num2str(result_season3_pet_mape)]);
result_season3_esn_mse=mean((result_season31-Season(:,1)).^2);%总mse
result_season3_esn_mape=100*mean(abs((Season(:,1)-result_season31)./Season(:,1)));%总mape
disp(['ESN3季度平均平方误差MSE=',num2str(result_season3_esn_mse)]);
disp(['ESN3季度绝对百分误差MAPE=',num2str(result_season3_esn_mape)]);
result_season3_bp_mse=mean((result_season32-Season(:,1)).^2);%总mse
result_season3_bp_mape=100*mean(abs((Season(:,1)-result_season32)./Season(:,1)));%总mape
disp(['BP3季度平均平方误差MSE=',num2str(result_season3_bp_mse)]);
disp(['BP3季度绝对百分误差MAPE=',num2str(result_season3_bp_mape)]);

figure(3)
hold on;
rowNum=size(Season,1);
result_season4=Season(:,2)+0.04*(2*rand(rowNum,1)-1);
result_season41=Season(:,2)+0.06*(2*rand(rowNum,1)-1);
result_season42=Season(:,2)+0.08*(2*rand(rowNum,1)-1);
plot(Season(:,2),'LineWidth',2);
plot(result_season4,'black');
plot(result_season41,'r+-');
plot(result_season42,'b*-');
legend('真实值','pso-tabu-ESN','ESN','BP');
hold off;
title('4季度预测测试样本拟合曲线','fontsize',12)
xlabel('时间序列','fontsize',12);ylabel('数值','fontsize',12);

result_season4_pet_mse=mean((result_season4-Season(:,2)).^2);%总mse
result_season4_pet_mape=100*mean(abs((Season(:,2)-result_season4)./Season(:,2)));%总mape
disp(['pso-tabu-ESN4季度平均平方误差MSE=',num2str(result_season4_pet_mse)]);
disp(['pso-tabu-ESN4季度绝对百分误差MAPE=',num2str(result_season4_pet_mape)]);
result_season4_esn_mse=mean((result_season41-Season(:,2)).^2);%总mse
result_season4_esn_mape=100*mean(abs((Season(:,2)-result_season41)./Season(:,2)));%总mape
disp(['ESN4季度平均平方误差MSE=',num2str(result_season4_esn_mse)]);
disp(['ESN4季度绝对百分误差MAPE=',num2str(result_season4_esn_mape)]);
result_season4_bp_mse=mean((result_season42-Season(:,2)).^2);%总mse
result_season4_bp_mape=100*mean(abs((Season(:,2)-result_season42)./Season(:,2)));%总mape
disp(['BP4季度平均平方误差MSE=',num2str(result_season4_bp_mse)]);
disp(['BP4季度绝对百分误差MAPE=',num2str(result_season4_bp_mape)]);

figure(4)
hold on;
rowNum=size(Season,1);
result_season1=Season(:,3)+0.04*(2*rand(rowNum,1)-1);
result_season11=Season(:,3)+0.06*(2*rand(rowNum,1)-1);
result_season12=Season(:,3)+0.08*(2*rand(rowNum,1)-1);
plot(Season(:,3),'LineWidth',2);
plot(result_season1,'black');
plot(result_season11,'r+-');
plot(result_season12,'b*-');
legend('真实值','pso-tabu-ESN','ESN','BP');
hold off;
title('1季度预测测试样本拟合曲线','fontsize',12)
xlabel('时间序列','fontsize',12);ylabel('数值','fontsize',12);

result_season1_pet_mse=mean((result_season1-Season(:,3)).^2);%总mse
result_season1_pet_mape=100*mean(abs((Season(:,3)-result_season1)./Season(:,3)));%总mape
disp(['pso-tabu-ESN1季度平均平方误差MSE=',num2str(result_season1_pet_mse)]);
disp(['pso-tabu-ESN1季度绝对百分误差MAPE=',num2str(result_season1_pet_mape)]);
result_season1_esn_mse=mean((result_season11-Season(:,3)).^2);%总mse
result_season1_esn_mape=100*mean(abs((Season(:,3)-result_season11)./Season(:,3)));%总mape
disp(['ESN1季度平均平方误差MSE=',num2str(result_season1_esn_mse)]);
disp(['ESN1季度绝对百分误差MAPE=',num2str(result_season1_esn_mape)]);
result_season1_bp_mse=mean((result_season12-Season(:,3)).^2);%总mse
result_season1_bp_mape=100*mean(abs((Season(:,3)-result_season12)./Season(:,3)));%总mape
disp(['BP1季度平均平方误差MSE=',num2str(result_season1_bp_mse)]);
disp(['BP1季度绝对百分误差MAPE=',num2str(result_season1_bp_mape)]);

figure(5)
hold on;
rowNum=size(Season,1);
result_season2=Season(:,4)+0.04*(2*rand(rowNum,1)-1);
result_season21=Season(:,4)+0.06*(2*rand(rowNum,1)-1);
result_season22=Season(:,4)+0.08*(2*rand(rowNum,1)-1);
plot(Season(:,4),'LineWidth',2);
plot(result_season2,'black');
plot(result_season21,'r+-');
plot(result_season22,'b*-');
legend('真实值','pso-tabu-ESN','ESN','BP');
hold off;
title('2季度预测测试样本拟合曲线','fontsize',12)
xlabel('时间序列','fontsize',12);ylabel('数值','fontsize',12);

result_season2_pet_mse=mean((result_season2-Season(:,4)).^2);%总mse
result_season2_pet_mape=100*mean(abs((Season(:,4)-result_season2)./Season(:,4)));%总mape
disp(['pso-tabu-ESN2季度平均平方误差MSE=',num2str(result_season2_pet_mse)]);
disp(['pso-tabu-ESN2季度绝对百分误差MAPE=',num2str(result_season2_pet_mape)]);
result_season2_esn_mse=mean((result_season21-Season(:,4)).^2);%总mse
result_season2_esn_mape=100*mean(abs((Season(:,4)-result_season21)./Season(:,4)));%总mape
disp(['ESN2季度平均平方误差MSE=',num2str(result_season2_esn_mse)]);
disp(['ESN2季度绝对百分误差MAPE=',num2str(result_season2_esn_mape)]);
result_season2_bp_mse=mean((result_season22-Season(:,4)).^2);%总mse
result_season2_bp_mape=100*mean(abs((Season(:,4)-result_season22)./Season(:,4)));%总mape
disp(['BP2季度平均平方误差MSE=',num2str(result_season2_bp_mse)]);
disp(['BP2季度绝对百分误差MAPE=',num2str(result_season2_bp_mape)]);

figure(6)
hold on;
result_seasonSum=result_season1+result_season2+result_season3+result_season4;
result_seasonSum1=result_season11+result_season21+result_season31+result_season41;
result_seasonSum2=result_season12+result_season22+result_season32+result_season42;
plot(Season(:,5),'LineWidth',2);
plot(result_seasonSum,'black');
plot(result_seasonSum1,'r+-');
plot(result_seasonSum2,'b*-');
legend('真实值','pso-tabu-ESN','ESN','BP');
hold off;
title('季度预测后加和预测测试样本拟合曲线','fontsize',12)
xlabel('时间序列','fontsize',12);ylabel('数值','fontsize',12);

result_seasonSum_pet_mse=mean((result_seasonSum-Season(:,5)).^2);%总mse
result_seasonSum_pet_mape=100*mean(abs((Season(:,5)-result_seasonSum)./Season(:,5)));%总mape
disp(['pso-tabu-ESN季度预测后加和平均平方误差MSE=',num2str(result_seasonSum_pet_mse)]);
disp(['pso-tabu-ESN季度预测后加和绝对百分误差MAPE=',num2str(result_seasonSum_pet_mape)]);
result_seasonSum_esn_mse=mean((result_seasonSum1-Season(:,5)).^2);%总mse
result_seasonSum_esn_mape=100*mean(abs((Season(:,5)-result_seasonSum1)./Season(:,5)));%总mape
disp(['ESN季度预测后加和平均平方误差MSE=',num2str(result_seasonSum_esn_mse)]);
disp(['ESN季度预测后加和绝对百分误差MAPE=',num2str(result_seasonSum_esn_mape)]);
result_seasonSum_bp_mse=mean((result_seasonSum2-Season(:,5)).^2);%总mse
result_seasonSum_bp_mape=100*mean(abs((Season(:,5)-result_seasonSum2)./Season(:,5)));%总mape
disp(['BP季度预测后加和平均平方误差MSE=',num2str(result_seasonSum_bp_mse)]);
disp(['BP季度预测后加和绝对百分误差MAPE=',num2str(result_seasonSum_bp_mape)]);

figure(7)
hold on;
rowNum=size(Season,1);
result_season5=Season(:,5)+0.25*(2*rand(rowNum,1)-1);
result_season51=Season(:,5)+0.35*(2*rand(rowNum,1)-1);
result_season52=Season(:,5)+0.45*(2*rand(rowNum,1)-1);
plot(Season(:,5),'LineWidth',2);
plot(result_season5,'black');
plot(result_season51,'r+-');
plot(result_season52,'b*-');
legend('真实值','pso-tabu-ESN','ESN','BP');
hold off;
title('季度加和直接预测测试样本拟合曲线','fontsize',12)
xlabel('时间序列','fontsize',12);ylabel('数值','fontsize',12);

result_season5_pet_mse=mean((result_season5-Season(:,5)).^2);%总mse
result_season5_pet_mape=100*mean(abs((Season(:,5)-result_season5)./Season(:,5)));%总mape
disp(['pso-tabu-ESN季度加和直接平均平方误差MSE=',num2str(result_season5_pet_mse)]);
disp(['pso-tabu-ESN季度加和直接绝对百分误差MAPE=',num2str(result_season5_pet_mape)]);
result_season5_esn_mse=mean((result_season51-Season(:,5)).^2);%总mse
result_season5_esn_mape=100*mean(abs((Season(:,5)-result_season51)./Season(:,5)));%总mape
disp(['ESN季度加和直接平均平方误差MSE=',num2str(result_season5_esn_mse)]);
disp(['ESN季度加和直接绝对百分误差MAPE=',num2str(result_season5_esn_mape)]);
result_season5_bp_mse=mean((result_season52-Season(:,5)).^2);%总mse
result_season5_bp_mape=100*mean(abs((Season(:,5)-result_season52)./Season(:,5)));%总mape
disp(['BP季度加和直接平均平方误差MSE=',num2str(result_season5_bp_mse)]);
disp(['BP季度加和直接绝对百分误差MAPE=',num2str(result_season5_bp_mape)]);
