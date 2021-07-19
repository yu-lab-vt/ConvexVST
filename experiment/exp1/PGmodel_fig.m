clc;clear;close all;
tic;
sigma = [5 10 15];
color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560];
for s = 1:3
dis_p = zeros(256,10001);
dis_g = normpdf(-500:500,0,sqrt(20));
dis = zeros(256,11001);
for ii = 0:255
    ii
    dis_p(ii+1,5001-sigma(s)*ii:10001) = poisspdf(0:sigma(s)*ii+5000,sigma(s)*ii);
    dis(ii+1,:) = conv(dis_p(ii+1,:),dis_g);
end
pdf = zeros(255,256);
for ii = 0:255
    % ii -> 3001
    pdf(ii+1,1) = sum(dis(ii+1,1:5501-ii));
    pdf(ii+1,256) = sum(dis(ii+1,5501-ii+255:end));
    pdf(ii+1,2:255) = dis(ii+1,5501-ii+1:5501-ii+254);
end
cdf = cumsum(pdf,2);
mean_ = zeros(1,256);
median_ = zeros(1,256);
mean_var = zeros(1,256);
median_var = zeros(1,256);
real_var = zeros(1,256);
for ii = 0:255
    mean_(ii+1) = sum(pdf(ii+1,:).*[0:255]);
    median_(ii+1) = find(cdf(ii+1,:)>0.5,1,'first') - 1;
    mean_var(ii+1) = sum(pdf(ii+1,:).*([0:255]-mean_(ii+1)).^2);
    median_var(ii+1) = sum(pdf(ii+1,:).*([0:255]-median_(ii+1)).^2);
    real_var(ii+1) = sum(pdf(ii+1,:).*([0:255]-ii).^2);
end
p(s) = plot(median_var,'Color', color(s,:) , 'LineWidth',2);hold on
plot(mean_var, 'Color', color(s,:), 'LineStyle', '--','LineWidth',2);
save(['dis1_' num2str(s) '.mat'],'pdf');
end
legend([p(1) p(2) p(3)],{'\chi = 5','\chi = 10','\chi = 15'},'FontSize',15);
xlabel('$\theta$','FontWeight','bold','FontSize',15,'Interpreter','latex');ylabel('$\sigma^2_{\theta}$','FontWeight','bold','FontSize',15,'Interpreter','latex');