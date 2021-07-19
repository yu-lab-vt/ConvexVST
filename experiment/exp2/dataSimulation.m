clc;clear;%close all
load data.mat
dataSize = size(data);
numFrame = dataSize(4);
realSignal = zeros(dataSize);
for ii = 1:numFrame
    realSignal(:,:,:,ii) = medfilt3(data(:,:,:,ii),[3 3 1]);     % most time consuming part
end

% % generate poisson-gaussian pdf
% % assume var = 8*x + 20 
cluster = discretize(realSignal,-0.5:255.5);
index = label2idx(cluster); 
sigma_p = 6;
sigma_g = sqrt(20);
dis = zeros(256,256);
variance = zeros(1,256);
data_simulation = zeros(dataSize);
for ii = 1:256
    ii
    len = length(index{ii});
    intensity = ii - 1;
    
    poiss_noise = poissrnd(sigma_p*intensity,len,1) - sigma_p*intensity;
    gauiss_noise = round(normrnd(0,sigma_g,[len 1]));
    
    data_temp = intensity + poiss_noise + gauiss_noise;
    data_temp(data_temp<0) = 0;
    data_temp(data_temp>255) = 255;
    distribution = histogram(data_temp,-0.5:255.5);
    dis(ii,:) = distribution.Values/len;
    variance(ii) = sum(dis(ii,:).*(([0:255]-intensity).^2));
    
    data_simulation(index{ii}) = data_temp;
end
save('simulation_6.mat','realSignal','data_simulation','dis','variance','-v7.3');
% for tt = 1:numFrame
%     ind = num2str(1000+tt); 
%     ind = ['.\dataSimulation\' ind(2:4)];
%     tifwrite(uint8(data_simulation(:,:,:,tt)),ind);
% end

