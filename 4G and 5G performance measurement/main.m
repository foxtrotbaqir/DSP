%% Created by Baqir Kazmi
% open to use for research and study purpose
% Load the CSV file
data_fl0 = readtable('fl2.csv');
% Extract relevant columns
timestamps = data_fl0.Time;
latitudes = data_fl0.Lat;
longitudes = data_fl0.Long;
throughput_5G_DL = data_fl0.NetPDSCHThp____PCell_;  % Net PDSCH Throughput (PDSCH Thp)
throughput_5G_UL = data_fl0.NetPUSCHThp____PCell_;  % Net PUSCH Throughput (PUSCH Thp)
throughput_4G_DL = data_fl0.PDSCHThrpt____PCell_;  % PDSCH Throughput (PDSCH Thp)
throughput_4G_UL = data_fl0.PUSCHThrpt____PCell_;  % PUSCH Throughput (PUSCH Thp)
rx_signal_power = data_fl0.RSRP____PCell_; % RSRP of Rx signal
rx_signal_quality = data_fl0.RSRQ____PCell_; % RSRQ of Rx signal 
rx_signal_SINR = data_fl0.SINRRx_0_____PCell_; % SINR of Rx signal
% rx_signal_BLER = data.BLER % yet to be decided
%% Create a figure
figure;

% Plot 5G NSA Downlink Throughput
subplot(2, 1, 1);
geoscatter(latitudes, longitudes, 10, throughput_5G_DL);  % Scatter plot with color-coded 5G DL throughput
colormap('jet');  % Choose colormap
c = colorbar;
c.Label.String = '5G NSA Net PDSCH Throughput (Mbps)';  % Set colorbar label
title('5G NSA Net PDSCH Throughput Heatmap');  % Set title
% Add annotations for axes
annotation('textbox',[0.5,0.05,0.1,0.05],'String','Longitude','EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.1,0.5,0.1,0.05],'String','Latitude','EdgeColor','none','HorizontalAlignment','center');

% Plot 4G Downlink Throughput
subplot(2, 1, 2);
geoscatter(latitudes, longitudes, 10, throughput_4G_DL);  % Scatter plot with color-coded 4G DL throughput
colormap('jet');  % Choose colormap
c = colorbar;
c.Label.String = '4G PDSCH Throughput (Mbps)';  % Set colorbar label
title('4G PDSCH Throughput Heatmap');  % Set title
% Add annotations for axes
annotation('textbox',[0.5,0.05,0.1,0.05],'String','Longitude','EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.1,0.5,0.1,0.05],'String','Latitude','EdgeColor','none','HorizontalAlignment','center');


%% % Create a geographical heatmap for 4G PUSCH Throughput
figure;
subplot(2, 1, 1);
geoscatter(latitudes, longitudes, 10, throughput_4G_UL);  % Scatter plot with color-coded PUSCH throughput
colormap('jet');  % Choose colormap
c = colorbar;
c.Label.String = '4G PUSCH Throughput (Mbps)';  % Set colorbar label
title('4G PUSCH Throughput Heatmap');  % Set title
annotation('textbox',[0.5,0.05,0.1,0.05],'String','Latitude','EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.1,0.5,0.1,0.05],'String','Longitude','EdgeColor','none','HorizontalAlignment','center');
% Create a 5G Net PUSCH Thorughput Heatmap
subplot(2, 1, 2);
geoscatter(latitudes, longitudes, 10, throughput_5G_UL);  % Scatter plot with color-coded PUSCH throughput
colormap('jet');  % Choose colormap
c = colorbar;
c.Label.String = '5G Net PUSCH Throughput (Mbps)';  % Set colorbar label
title('5G Net PUSCH Throughput Heatmap');  % Set title
annotation('textbox',[0.5,0.05,0.1,0.05],'String','Latitude','EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.1,0.5,0.1,0.05],'String','Longitude','EdgeColor','none','HorizontalAlignment','center');
%% RSRP and RSRQ
%% Create Geographical Heatmap for RSRP
figure;
subplot(2, 1, 1);
geoscatter(latitudes, longitudes, 10, rx_signal_power);  % Scatter plot with color-coded RSRP values
colormap('jet');  % Choose colormap
c = colorbar;
c.Label.String = 'RSRP (dBm)';  % Set colorbar label
title('RSRP Heatmap');  % Set title
annotation('textbox',[0.5,0.05,0.1,0.05],'String','Latitude','EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.1,0.5,0.1,0.05],'String','Longitude','EdgeColor','none','HorizontalAlignment','center');

% Create Geographical Heatmap for RSRQ
subplot(2, 1, 2);
geoscatter(latitudes, longitudes, 10, rx_signal_quality);  % Scatter plot with color-coded RSRQ values
colormap('jet');  % Choose colormap
c = colorbar;
c.Label.String = 'RSRQ';  % Set colorbar label
title('RSRQ Heatmap');  % Set title
annotation('textbox',[0.5,0.05,0.1,0.05],'String','Latitude','EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.1,0.5,0.1,0.05],'String','Longitude','EdgeColor','none','HorizontalAlignment','center');
%% Create Histogram for SINR, RSRQ, and DL Throughput (5G and 4G)
figure;
subplot(4, 1, 1);
histogram(rx_signal_SINR);
xlabel('SINR');
ylabel('Frequency');
title('Histogram of SINR');

subplot(4, 1, 2);
histogram(rx_signal_quality);
xlabel('RSRQ');
ylabel('Frequency');
title('Histogram of RSRQ');

subplot(4, 1, 3);
histogram(throughput_5G_DL);
xlabel('5G DL (Mbps)');
ylabel('Frequency');
title('Histogram of 5G DL Throughput');

subplot(4,1,4);
histogram(throughput_4G_DL);
xlabel('4G DL (Mbps)');
ylabel('Frequency');
title('Histogram of 4G DL Throughput');

