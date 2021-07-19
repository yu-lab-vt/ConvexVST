clc;clear;close all;
dbstop if error
addpath('.\bm3d');
addpath('C:\Program Files\Mosek\9.2\toolbox\R2015a');
dataset_path = '..\exp3\dataset';
dataset_list = dir(dataset_path);
result_error = zeros(12,4);
result_bar = zeros(12,4);
PSNR_all = zeros(4,20,12);
MSE_all = zeros(4,20,12);
for dd = 1:12
    close all;
dataset = dataset_list(dd+2).name;
PSNR_no = zeros(1,20);
PSNR_ori = zeros(1,20);  MSE_ori = PSNR_ori;
PSNR_convex = zeros(1,20); MSE_convex = PSNR_convex;
PSNR_gat = zeros(1,20);  MSE_gat = PSNR_gat;
PSNR_foi = zeros(1,20);  MSE_foi = PSNR_foi;
% only use raw data
data = zeros(512,512,20,50); % resolution FOV sample
distribution = zeros(1,256);

for ii = 1:20
    filepath = [dataset_path '\' dataset '\raw\' num2str(ii)];
    name_list = dir(filepath);
    for jj = 1:50
        data(:,:,ii,jj) = imread([filepath '\' name_list(jj+2).name]);
    end
end

% autoclipping
clip_number = 1000;
histEdges = 0.5:254.5;
cluster = discretize(median(data,4),histEdges);
index = label2idx(cluster);
left_number = 0;
left_intensity = 1;
right_number = 0;
right_intensity = 254;
for ii = 1:254
    left_number = left_number + length(index{ii});
    if left_number >= clip_number
        left_intensity = ii;
        break
    end
end
for ii = length(index):-1:1
    right_number = right_number + length(index{ii});
    if right_number >= clip_number
        right_intensity = ii;
        break
    end
end
data(data<left_intensity) = left_intensity;
data(data>right_intensity) = right_intensity;
data = data-left_intensity; 
right_intensity = right_intensity - left_intensity;
histEdges = 0.5:right_intensity-0.5;
binEdges = -0.5:right_intensity+0.5;
histCenters = (histEdges(1:end-1) + histEdges(2:end))/2;
binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;
realSignal = median(data(:,:,:,1:25),4);
realSignal = repmat(realSignal,1,1,1,25);


% our method
figure(2);
[variance,his,parameters] = histogramCount(data(:,:,:,1:25),realSignal,binEdges,binEdges);   % generate variance curve
realSignal = median(data(:,:,:,26:50),4);
realSignal = repmat(realSignal,1,1,1,25);
variance_original = variance;
parameters_original = parameters;
[PSNR_ori, MSE_ori, data_denoised] = PSNR_check2(data(:,:,:,26:50), realSignal, variance, binCenters);
figure(1);plot(parameters.histCenters+left_intensity,variance,'o');hold on;
figure(2);
[variance,his,parameters] = histogramCount(data(:,:,:,1:25),realSignal,histEdges,binEdges);  % generate histogram wihout ends
[variance, stabilizeFunction] = optimization(his,parameters);
figure(3);plot(0:right_intensity,stabilizeFunction,'LineWidth',2);hold on

sf_old = stabilizeFunction;
stabilizeFunction(1) = stabilizeFunction(2);
stabilizeFunction(end) = stabilizeFunction(end-1);
stabilizeFunction = (stabilizeFunction - min(stabilizeFunction))/(max(stabilizeFunction) - min(stabilizeFunction));
stabilizeFunction = stabilizeFunction*(max(sf_old) - min(sf_old)) + min(sf_old);

realSignal = median(data(:,:,:,26:50),4);
realSignal = repmat(realSignal,1,1,1,25);
realSignal = realSignal+1;
realSignal = stabilizeFunction(realSignal);
data_stab = data+1;
data_stab = stabilizeFunction(data_stab);
figure(2);
[variance,his,parameters] = histogramCount(data_stab(:,:,:,26:50),realSignal,binEdges,binEdges);  % validation
realSignal = median(data(:,:,:,26:50),4);
realSignal = repmat(realSignal,1,1,1,25);
[PSNR_convex, MSE_convex, data_denoised] = PSNR_check(data(:,:,:,26:50), realSignal, sqrt(median(variance)), binCenters, stabilizeFunction);

% generalized anscome transformation
start_ind = floor(length(variance_original)*0.1);
end_ind = ceil(length(variance_original)*0.5);
linear_fit = fit(parameters_original.histCenters(start_ind:end_ind)',variance_original(start_ind:end_ind)','poly1');
poisson = sqrt(max(linear_fit.p1,0)); gaussian = sqrt(max(linear_fit.p2,0));
stabilizeFunction = 2/poisson*sqrt(poisson*parameters_original.binCenters+3*poisson^2/8+gaussian^2);
c = (max(stabilizeFunction) - min(stabilizeFunction));
stabilizeFunction = (stabilizeFunction - min(stabilizeFunction))/(max(stabilizeFunction) - min(stabilizeFunction))...
                        * (max(parameters_original.binCenters) - min(parameters_original.binCenters));
figure(3);plot(0:right_intensity,stabilizeFunction,'LineWidth',2);              
realSignal = median(data(:,:,:,26:50),4);
realSignal = repmat(realSignal,1,1,1,25);
realSignal = realSignal+1;
realSignal = stabilizeFunction(realSignal);
data_stab = data+1;
data_stab = stabilizeFunction(data_stab);
figure(2);
[variance,his,parameters] = histogramCount(data_stab(:,:,:,26:50),realSignal,binEdges,binEdges);  % validation
realSignal = median(data(:,:,:,26:50),4);
realSignal = repmat(realSignal,1,1,1,25);
[PSNR_gat, MSE_gat, data_denoised] = PSNR_check(data(:,:,:,26:50), realSignal, c, binCenters, stabilizeFunction);

% foi
realSignal = median(data(:,:,:,1:25),4);
realSignal = repmat(realSignal,1,1,1,25);
figure(2);
[variance,his,parameters] = histogramCount(data(:,:,:,1:25),realSignal,binEdges,binEdges);   % generate variance curve
ou = 1.5; ol = 1.5; ru = 0.8; rl = 0.8;
ru1 = 0.2; rl1 = 0.2; ru2 = 0.5; rl2 = 0.5;
sigmaf = sqrt(variance_original);
fk = 0:right_intensity;
for kk = 1:500
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
fk = fk/max(fk)*right_intensity;
for ii = 1:right_intensity+1
    variance(ii) = sum(his(ii,:).*(fk- fk(ii)).^2);
end
sigmaf = sqrt(variance);
end
figure(3);plot(0:right_intensity,fk);
realSignal = median(data(:,:,:,26:50),4);
realSignal = repmat(realSignal,1,1,1,25);
realSignal = realSignal+1;
realSignal = fk(realSignal);
data_stab = data+1;
data_stab = fk(data_stab);
figure(2);
[variance,his,parameters] = histogramCount(data_stab(:,:,:,26:50),realSignal,binEdges,binEdges);  % validation
realSignal = median(data(:,:,:,26:50),4);
realSignal = repmat(realSignal,1,1,1,25);
[PSNR_foi, MSE_foi, data_denoised] = PSNR_check(data(:,:,:,26:50), realSignal, c, binCenters, fk);
% writeDenoised(['D:\dropbox\VSTApplication\Exp-Comparison\exp2\denoising\' dataset '\foi'], data_denoised);
PSNR_all(:,:,dd) = [PSNR_ori; PSNR_convex; PSNR_gat; PSNR_foi;];
MSE_all(:,:,dd) = [MSE_ori; MSE_convex; MSE_gat; MSE_foi;];
end
save('PSNR','PSNR_all', 'MSE_all');

function [PSNR, MSE, data_denoised] = PSNR_check(data, realSignal, sigma, binCenters, stabilizeFunction)
I_MAX = 255;
% stabilization
data = interp1(binCenters,stabilizeFunction,data,'linear','extrap');
data = data/I_MAX;
sigma = sigma/I_MAX;
disp('Denoising started')
PSNR = zeros(1,20);
MSE = zeros(1,20);
[x, y, t, z] = size(data);
tic;
data_denoised = zeros(size(data));
for tt = 1:20
    tt
    data_est = zeros(x,y,1,z);
    parfor zz = 1:z
        data_est(:,:,1,zz) = BM3D(data(:,:,tt,zz), sigma);
    end
    data_est = data_est*I_MAX;
    data_est = interp1(stabilizeFunction(2:end-1), binCenters(2:end-1),data_est,'linear','extrap');
    data_est(data_est < stabilizeFunction(2)) = 0;
    data_est(data_est > stabilizeFunction(end-1)) = stabilizeFunction(end-1);
    MSE(tt) = mean(((realSignal(:,:,tt,:)-data_est(:,:,1,:)).^2).*(realSignal(:,:,tt,:)>5),'all');
    PSNR(tt) = 10*log10(I_MAX^2/MSE(tt));
    data_denoised(:,:,tt,:) = data_est(:,:,1,:);
end
time = toc;
fprintf('Denoising completed %.1fs: PSNR %.2fdB \n', time, PSNR);
end

function [PSNR, MSE, data_denoised] = PSNR_check2(data, realSignal, variance, binCenters)
I_MAX = 255;
[x, y, t, z] = size(data);
sigma = zeros(size(data));
for tt = 1:20
    for zz = 1:z
        sigma(:,:,tt,zz) = medfilt2(data(:,:,tt,zz),[3 3]);
    end
end
sigma = interp1(binCenters, variance, sigma);
data = data/I_MAX;
disp('Denoising started')
PSNR = zeros(1,20);
MSE = zeros(1,20);
data_denoised = zeros(size(data));

tic;
for tt = 1:20
    tt
    data_est = zeros(x,y,1,z);
    parfor zz = 1:z
        data_est(:,:,1,zz) = BM3D(data(:,:,tt,zz), sigma(:,:,tt,zz));
    end
    data_est = data_est*I_MAX;
    MSE(tt) = mean(((realSignal(:,:,tt,:)-data_est(:,:,1,:)).^2).*(realSignal(:,:,tt,:)>5),'all');
    PSNR(tt) = 10*log10(I_MAX^2/MSE(tt));
    data_denoised(:,:,tt,:) = data_est(:,:,1,:);
end
time = toc;
fprintf('Denoising completed %.1fs: PSNR %.2fdB \n', time, PSNR);
end


function [variance, stabilizeFunction] = optimization(his,parameters)
if any(sum(his,2) < 1 - 1e-6) || any(sum(his,2) > 1 + 1e-6)
    warning('The sum of some histograms may be not 1.');
    his = his./sum(his,2);
end
if ~isfield(parameters,'histCenters')
error('parameters.histCenters is necessary.');
end
if ~isfield(parameters,'binCenters')
error('parameters.binCenters is necessary.');
end
histCenters = parameters.histCenters;
binCenters = parameters.binCenters; 
if min(histCenters) < min(binCenters) || max(histCenters) > max(binCenters)
    error('The range of histCenters should not exceed the range of binCenters.');
end
numOfHists = length(histCenters);
numOfBins = length(binCenters);
if ~isfield(parameters,'var_ratio')
    var_ratio = ones(1,numOfHists);
else
    var_ratio = parameters.var_ratio;
end
if ~isfield(parameters,'step')
    step = 1;
else
    step = parameters.step;
end
if ~isfield(parameters,'algorithm')
    algorithm = 'aaaaa';
else
    algorithm = 'bbbbb';
end

% optimization part first step
tic;
disp('Start First Step...');

% MOSEK Part`
Hess = cell(1,numOfHists);
for ii = 1:numOfHists
    Hess{ii} = zeros(numOfBins-1,numOfBins-1);
    left_ind = find(binCenters<histCenters(ii),1,'last'); % also the split variable
    right_ind = left_ind+1; 
    alpha = (histCenters(ii) - binCenters(left_ind))/(binCenters(right_ind) - binCenters(left_ind));
    beta = (binCenters(right_ind) - histCenters(ii))/(binCenters(right_ind) - binCenters(left_ind));
    if isempty(left_ind)
        left_ind = 0; right_ind = 1; alpha = 1; beta = 0;
    end
    for jj = 1:numOfBins
        if jj < left_ind
            Hess{ii}(jj:left_ind-1,jj:left_ind-1) = Hess{ii}(jj:left_ind-1,jj:left_ind-1) + his(ii,jj);
            Hess{ii}(left_ind,jj:left_ind-1) = Hess{ii}(left_ind,jj:left_ind-1) + alpha*his(ii,jj);
            Hess{ii}(jj:left_ind-1,left_ind) = Hess{ii}(jj:left_ind-1,left_ind) + alpha*his(ii,jj);
            Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + alpha^2*his(ii,jj);
        elseif jj == left_ind
            Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + alpha^2*his(ii,jj);
        elseif jj == right_ind
            if left_ind > 0
                Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + beta^2*his(ii,jj);
            end
        elseif jj > right_ind
            Hess{ii}(right_ind:jj-1,right_ind:jj-1) = Hess{ii}(right_ind:jj-1,right_ind:jj-1) + his(ii,jj);
            if left_ind > 0 
                Hess{ii}(left_ind,right_ind:jj-1) = Hess{ii}(left_ind,right_ind:jj-1) + beta*his(ii,jj);
                Hess{ii}(right_ind:jj-1,left_ind) = Hess{ii}(right_ind:jj-1,left_ind) + beta*his(ii,jj);
                Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + beta^2*his(ii,jj);
            end
        end
    end
    Hess{ii} = tril(Hess{ii})*2/var_ratio(ii);
    Hess{ii}(Hess{ii}<1e-16) = 0;
end
prob.qcsubk = [];
prob.qcsubi = [];
prob.qcsubj = [];
prob.qcval = [];
for ii = 1:numOfHists
    [row,column] = find(Hess{ii});
    prob.qcsubk = [prob.qcsubk; ii*ones(length(row),1)];
    prob.qcsubi = [prob.qcsubi; row];
    prob.qcsubj = [prob.qcsubj; column];
    prob.qcval = [prob.qcval; Hess{ii}(sub2ind(size(Hess{ii}),row,column))];
end
prob.c = [zeros(1,numOfBins-1) 1];
prob.buc = [zeros(numOfHists,1); binCenters(end)-binCenters(1)];
prob.blc = [-inf*ones(numOfHists,1); binCenters(end)-binCenters(1)];
prob.blx = zeros(numOfBins,1);
prob.a = repmat([zeros(1,numOfBins-1) -1],[numOfHists 1]);
prob.a = [prob.a; [ones(1,numOfBins-1) 0]];
param.MSK_IPAR_INTPNT_SOLVE_FORM = 'MSK_SOLVE_DUAL';
[~,res] = mosekopt('minimize',prob,param);

x = res.sol.itr.xx(1:end-1)';
stabilizeFunction = cumsum([binCenters(1) x]);
f = @(x)objectiveFunction(x,his,numOfHists,binCenters,histCenters);
variance = f(x);
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

function inv_fx = inverseFunction(fx)

% For any monotonic increasing function f(x) in [0-255], this function will output
% f-1(x)

len = length(fx);

inv_fx = zeros(1,len);

for ii = 2:len
    if round(fx(ii)) >= ceil(fx(ii-1))
        for jj = ceil(fx(ii-1)):round(fx(ii))
            if jj > 0
                inv_fx(jj+1) = (jj - fx(ii-1))/(fx(ii) - fx(ii-1)) + ii - 2;
            end
        end
    end
end
end