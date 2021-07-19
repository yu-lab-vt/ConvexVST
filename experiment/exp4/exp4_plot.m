clc;clear;close all;
load PSNR.mat
dataset_path = 'D:\dropbox\VSTApplication\Exp-Comparison\exp2\dataset';
dataset_list = dir(dataset_path);
pp = zeros(1,12);
for dd = 1:12
    temp = mean(PSNR_all(:,:,dd)');
    pp(dd) = (temp(2) - max(temp(3), temp(4)))/(max(temp) - min(temp));
end
[a, pp] = sort(-pp);
for dd = 1:12
    subplot(3,4,dd);
    p = boxplot(PSNR_all(:,:,pp(dd))','colors','krmb','Whisker',10);
    title(dataset_list(pp(dd)+2).name,'Interpreter','none','FontName','Times New Roman','FontSize',12);
    set(gca,'XTick',[]);            
    hold on
end
axis([0.5 4.5 25 35])  
hLegend = legend(flip(findall(gca,'Tag','Box')), {'No VST','ConvexOpt','GAT','Foi'});
10*log10(255^2/mean(MSE_all(1,:,:),'all'))