clc;clear;close all;
% dbstop if error

result_error1 = zeros(3,3);
result_bar = zeros(3,3);
for dd = 1:3
close all;

load(['dis1_' num2str(dd) '.mat']);
% load('result1_1.mat','variance');

parameters.binCenters = 0:255;
parameters.histCenters = 1:254;
variance = zeros(1,254);
for ii = 1:256
    variance(ii) = sum(pdf(ii,:).*([0:255]- ii-1).^2);
end
variance_original = variance;
parameters_original = parameters;
load(['result1_' num2str(dd) '.mat']);
figure(1);plot(parameters.binCenters,[0 cumsum(x)],'LineWidth',2);hold on;
% [variance, stabilizeFunction] = optimization(pdf(2:255,:),parameters);
figure(2);plot(parameters.binCenters,variance,'LineWidth',2);hold on;
result_error1(dd,1) = result_error(6)/256;
result_bar(dd,1) = max(variance) - min(variance);

% generalized GAT transformation
start_ind = floor(length(variance_original)*0.1);
end_ind = ceil(length(variance_original)*0.5);
linear_fit = fit(parameters_original.histCenters(start_ind:end_ind)',variance_original(start_ind:end_ind)','poly1');
poisson = sqrt(linear_fit.p1); gaussian = sqrt(max(linear_fit.p2,0));
stabilizeFunction = 2/poisson*sqrt(poisson*parameters_original.binCenters+3*poisson^2/8+gaussian^2);
c = (max(stabilizeFunction) - min(stabilizeFunction));
stabilizeFunction = (stabilizeFunction - min(stabilizeFunction))/(max(stabilizeFunction) - min(stabilizeFunction))...
                        * (max(parameters_original.binCenters) - min(parameters_original.binCenters));
for ii = 1:256
    variance(ii) = sum(pdf(ii,:).*(stabilizeFunction- stabilizeFunction(ii)).^2);
end
figure(1);plot(parameters.binCenters,stabilizeFunction,'LineWidth',2);hold on;
figure(2);plot(parameters.binCenters,variance,'LineWidth',2);hold on;
result_error1(dd,2) = sum(abs(variance-median(variance)))/256;
result_bar(dd,2) = max(variance) - min(variance);

% foi method
ou = 1.5; ol = 1.5; ru = 0.8; rl = 0.8;
ru1 = 0.2; rl1 = 0.2; ru2 = 0.5; rl2 = 0.5;

sigmaf = sqrt(variance_original);
fk = 0:255;

for kk = 1:150
ef = sigmaf - c;
ef_bar = max(-rl2,min(ru2,ef));
sigmaf_bar = max(c-rl2,min(c+ru2,sigmaf));
x_0lf = ef_bar<0;
x_0rg = ef_bar>0;
x_ru1lf = ef_bar<=ru1;
x_ru1rg = ef_bar>ru1;
x_rl1lf = ef_bar<=-rl1;
x_rl1rg = ef_bar>-rl1;

phi_ef = ru*x_0rg.*( (1-((ef-ru1)/ru1).^2).^(ou-1).*x_ru1lf + x_ru1rg )...
      +  rl*x_0lf.*( (1-((ef+rl1)/rl1).^2).^(ol-1).*x_rl1rg + x_rl1lf );
ifk = 1 - phi_ef.*ef_bar./sigmaf_bar;
increment = (ifk(1:end-1)+ifk(2:end))/2;
fk = [ 0 cumsum(increment.*diff(fk))];



for ii = 1:256
    variance(ii) = sum(pdf(ii,:).*(fk- fk(ii)).^2);
end
sigmaf = sqrt(variance);
end
fk = fk/max(fk)*255;
for ii = 1:256
    variance(ii) = sum(pdf(ii,:).*(fk- fk(ii)).^2);
end
figure(1); plot(fk,'LineWidth',2);hold on
plot(parameters.binCenters,0:255,'--','LineWidth',2);
figure(2); plot(variance,'LineWidth',2);hold on;
plot(variance_original,'--','LineWidth',2);
result_error1(dd,3) = sum(abs(variance-median(variance)))/256;
result_bar(dd,3) = max(variance) - min(variance);
figure(1);xlabel('\theta','FontWeight','bold','FontSize',15);
ylabel('$f(\theta)$','FontWeight','bold','Interpreter','latex','FontSize',15);
legend('ConvexOpt','GAT','Foi','Original','Location','northwest');

figure(2);xlabel('\theta','FontWeight','bold','FontSize',15);
ylabel('$\sigma^2_{\theta}$','FontWeight','bold','Interpreter','latex','FontSize',15);
legend('ConvexOpt','GAT','Foi','Original','Location','northwest');

end

function F = objectiveFunction(x,his,numOfHists,binCenters,histCenters)
    for ii = 1:numOfHists
        left_ind = find(binCenters<histCenters(ii),1,'last');
        if isempty(left_ind)
            F(ii) = sum(his(ii,:).*([0 cumsum(x)]).^2);
        else
            beta = (binCenters(left_ind+1) - histCenters(ii))/(binCenters(left_ind+1) - binCenters(left_ind));
            F(ii) = sum(his(ii,:).*([0 cumsum(x)]-(sum(x(1:left_ind))-beta*x(left_ind))).^2);
        end
    end
end

function [variance,his,parameters] = histogramCount(data,realSignal,histEdges,binEdges)
    data = data(:);                  
    realSignal = realSignal(:);
    
    cluster = discretize(realSignal,histEdges);
    index = label2idx(cluster); 

    numOfHists = length(histEdges) - 1;
    numOfBins = length(binEdges) - 1;
    histCenters = (histEdges(1:end-1) + histEdges(2:end))/2;
    binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;

    if min(cluster) > 1
        warning(['The smallest histogram does not contain any samples. Please'...
            ' double check options.histEdges.']);
    end
    if max(cluster) < length(histCenters)
        warning(['The largest histogram does not contain any samples. Please'...
            ' double check options.histEdges.']);
        histEdges = histEdges(1:length(index)+1);
        histCenters = histCenters(1:length(index));
    end

    % merge empty histograms
    index_label = zeros(1,length(index));
    for ii = 1:length(index)-1
        if length(index{ii}) < 200
            index_label(ii) = 1;
            if length(index{ii+1}) + length(index{ii}) > 0 
                histCenters(ii+1) = (histCenters(ii+1)*length(index{ii+1}) + ...
                    histCenters(ii)*length(index{ii}))/(length(index{ii+1})+length(index{ii}));
            end
            index{ii+1} = [index{ii+1};index{ii};];
        end
    end
    histEdges(logical([0 index_label])) = []; 
    histCenters(logical(index_label)) = [];
    index = index(logical(1-index_label));

    variance = zeros(1,length(index));
    his = zeros(length(index),numOfBins);
    for ii = 1:length(index)
        noise = data(index{ii}) - realSignal(index{ii});
        distribution = histogram(noise + histCenters(ii), binEdges);
        his(ii,:) = distribution.Values/length(noise);
        variance(ii) = sum(his(ii,:).*(binCenters-histCenters(ii)).^2);
    end
    parameters.histCenters = histCenters;
    parameters.histEdges = histEdges;
    parameters.binCenters = binCenters;
    parameters.binEdges = binEdges;
end
