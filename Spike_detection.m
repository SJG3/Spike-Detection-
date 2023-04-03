%% 
clear
clc 
%Sylvester James Gates III, Losert Lab, University of Maryland College Park
%2021-23
%This code is used to find peaks/spikes from a deltaF/F matrix and then can
%be used to output figures which include bar graphs for the active cells in
%population, percentage of active cells, and percetage of spikes in active
%cells.
%% conditions
tot_time_m = 3.333; %total time in minutes
fps = 5;            % FPS of video used

%% Use FindPeaks to automate detecting spikes, spike width at half prominence, and height

%create empty vectors which will store the spikes widths and prominence
widths_spikes = NaN(size(DFF,1));
prom_spikes = NaN(size(DFF,1));

%for every cell, use findpeaks with moving mean (10) to find peaks with
%prominence greater than 20. store the spikes widths/ppriminence of each cell on different
%row, with each spikes width/prominence on a different column
zen = DFF*0;

spike_width = 10; % wanted spike width in seconds
    %initally use 30seconds, then 10seconds, then 3seconds
for n = 1:size(DFF,1)
    [pks,locs,w,p] = findpeaks(movmean(DFF(n,:)',1),'MinPeakProminence',10,'MaxPeakWidth',fps*spike_width);
    number_of_spikes(n) = size(locs,1);
    widths_spikes(n,1:size(w',2)) = [w'];
    prom_spikes(n,1:size(w',2)) = [p'];
    
    for n_l = 1:size(locs,1)
        zen(n,locs(n_l)) = 1;
    end
end 

%use the stored widths and prominence to get a mean spike width for each
%cell and prominence
widths_spikes_mean_per_cell = nanmean(widths_spikes,2);
prom_spikes_mean_per_cell = nanmean(prom_spikes,2);

%call active cells those with greater than 3 spikes
active_cells_spikes = number_of_spikes(number_of_spikes>3); 
active_cells_mean_spikes = mean(active_cells_spikes);
%% Percentage cells active (ie > 3 spikes)
total_cells = size(DFF,1);
pct = size(active_cells_spikes,2) / total_cells
tot_spike_num = sum(number_of_spikes) 
tot_active_spike_num = sum(active_cells_spikes)
spikes_per_min = tot_spike_num/tot_time_m
%% Histogram of count as spike per min
figure(5347259);
histogram (number_of_spikes/tot_time_m , 'Normalization','probability', 'BinWidth',0.1, 'FaceAlpha',0.15);
hold on; 
%% FIGURES
figure(26543245);
subplot(1,2,1); scatter(prom_spikes_mean_per_cell,1:size(DFF,1),'.'); xlabel("Mean Spike Prominence (Intensity)"); ylabel("Cell Index"); hold on;
subplot(1,2,2); scatter(widths_spikes_mean_per_cell,1:size(DFF,1),'.');  hold on; 
% errorbar(widths_spikes_mean_per_cell,1:size(DFF,1),ans,'horizontal','.');
xlabel("Mean Spike Width (frames)"); ylabel("Cell Index");
%%      NOT NEEDED - TESTING for spike finder
figure(45432); 
nom = 35; 
 
for nom = 1:size(DFF,1)%[1:10:size(DFF,1)]
    % findpeaks(detrend(movmean(DFF(nom,:)',10),3),'MinPeakProminence',10,'Annotate','extents');title("cell "+nom);
    findpeaks(movmean(DFF(nom,:)',fps),'MinPeakProminence',15,'MinPeakWidth',fps*0.5,'MaxPeakWidth',fps*spike_width,'Annotate','extents');
%     findpeaks(movmean(DFF(nom,:)',fps),'MinPeakProminence', mean(reshape(abs(diff(DFF_smooth)),[],1))* 3 ,'MinPeakWidth',fps*0.5,'MaxPeakWidth',fps*spike_width,'Annotate','extents');
%     findpeaks(movmean(DFF(nom,:)',fps),'MinPeakProminence', std(reshape(DFF_smooth(:,:),[],1))*3,'MinPeakWidth',fps*0.5,'MaxPeakWidth',fps*spike_width,'Annotate','extents');
   
ylim([- max(DFF(:)) ,  max(DFF(:))]);
%     hold on
pause(0.2);
end 
%%      NOT NEEDED TESTING is detrending needed? 
nom = 52;
y = detrend(movmean(DFF(nom,:)',10),3);

figure(3453);
subplot(1,3,1); plot(y); ylim([- max(DFF(:)) ,  max(DFF(:))]);

subplot(1,3,2); plot(movmean(DFF(nom,:)',10));ylim([- max(DFF(:)) ,  max(DFF(:))]);

subplot(1,3,3); plot(smoothdata(DFF(nom,:)','sgolay',10));ylim([- max(DFF(:)) ,  max(DFF(:))]);
%% Use spike widths and prominence and scatter plot

% take prominence/widths and linearize the vector
lin_prom = prom_spikes';
lin_prom = lin_prom(:);

lin_width = widths_spikes';
lin_width = lin_width(:);

%remove NaNs from the vector
lin_prom = lin_prom(~isnan(lin_prom));
lin_width = lin_width(~isnan(lin_width));

%create cell labels which follow the linearized versions, to be used as
%colormap in plot
for n = 1:size(DFF,1)
    labss(n,:) = prom_spikes(n,:) * 0 + n;
end 
labss_d = labss';
labss_d = labss_d(:);
labss_d = labss_d(~isnan(labss_d));
%% FIGURES
%plot spikeprominence vs width
figure(4567);
scatter (lin_prom,lin_width,[],labss_d,'filled'); 
colormap(turbo); 
a  = colorbar;
a.Label.String = 'Cell Index';

xlabel("Prominence (Intensity)");
ylabel("Width (half prominence)");
%%      NOT NEEDED total number of cells/spikes vs prominence or 1/duration
% 
% subplot(1,2,1); plot(total_number_cells_p); xlabel("Prominence"); ylabel("Total Number of Cells ")
% subplot(1,2,2); plot(total_number_cells_w); xlabel("Width"); ylabel("Total Number of Cells")
%% total number of cells/spikes vs prominence or 1/duration

for val = 1:max(lin_prom(:))
    hol = lin_prom>=val;
    total_number_spikes_p(1,val) = sum(hol(:) == 1);
end 
    nrm_tns_p = total_number_spikes_p/max(total_number_spikes_p);

n = 99;
tot_n = n + 1;
refs = 0: max(1./lin_width(:))/99 :  max(1./lin_width(:));

for val = 1:tot_n
    hol = (1 ./ lin_width) >= refs(val);
    total_number_spikes_w(1,val) = sum(hol(:) == 1);
end 
    nrm_tns_w = total_number_spikes_w/max(total_number_spikes_w);
%% FIGURES
figure(2546);
subplot(2,2,1); plot(nrm_tns_p); xlabel("Prominence"); ylabel("% Total Number of Spikes > Prominence "); hold on; 
subplot(2,2,2); plot(refs,nrm_tns_w); xlabel("1/Width(frame duration)"); ylabel("% Total Number of Spikes > 1/Width(frame duration)"); hold on; 
    
subplot(2,2,3); plot(total_number_spikes_p); xlabel("Prominence"); ylabel("Total Number of Spikes > Prominence"); hold on; 
subplot(2,2,4); plot(refs,total_number_spikes_w); xlabel("Width"); ylabel("Total Number of Spikes"); hold on; 

%% BAR GRAPHS - percetnage active cells in conditions 

%for 30s max width
% pct_spikes_ctrl = [0 , 0.0482 , 0.0154];
% pct_spikes_AAm = [0.1357 , 0.2473 , 0.1147, 0.0343 ];
% pct_spikes_AATy = [0.5344 , 0.7116 , 0.6632];

%for 10s max width
% pct_spikes_ctrl = [0 , 0.0044, 0.0330];
% pct_spikes_AAm = [0.1181 , 0.1702 , 0.0800 , 0.0252 ];
% pct_spikes_AATy = [0.4904 , 0.6667 , 0.5953];

%for 3s max width
pct_spikes_ctrl = [0 , 0 , 0];
pct_spikes_AAm = [0 , 0 , 0 , 0.0046];
pct_spikes_AATy = [0.0579 , 0.1111 , 0.0574];

x = categorical({'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});
x = reordercats(x,{'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});
data = [mean(pct_spikes_ctrl),mean(pct_spikes_AAm),mean(pct_spikes_AATy)]';
errhigh = [std(pct_spikes_ctrl),std(pct_spikes_AAm),std(pct_spikes_AATy)];
errlow  = [std(pct_spikes_ctrl),std(pct_spikes_AAm),std(pct_spikes_AATy)];

figure(254); subplot(1,3,1);
b = bar(x,data);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.14 0.14 0.8];
    b.CData(2,:) = [0.8 0.14 0.14];
    b.CData(3,:) = [0.14 0.8 0.14];
ylabel("Percentage of Active Cells (>3 spikes)");
ylim([-0.1 1.1])

hold on
    er = errorbar(x,data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
hold off

%stats
[h,p] = ttest2(pct_spikes_ctrl , pct_spikes_AAm , 'Vartype' , 'unequal')
[h,p] = ttest2(pct_spikes_ctrl , pct_spikes_AATy , 'Vartype' , 'unequal')
[h,p] = ttest2(pct_spikes_AAm , pct_spikes_AATy , 'Vartype' , 'unequal')
%% BAR GRAPHS - spikes per min (Freuqency) in conditions 

%for 30s max width
% spikes_per_min_ctrl = [2.5 , 22.2222, 27.25];
% spikes_per_min_AAm = [74 , 93.5556 , 52 , 46];
% spikes_per_min_AATy = [318.2222 , 532.5556 , 425];

%for 10s max width
% spikes_per_min_ctrl = [0 , 12 , 15.8889];
% spikes_per_min_AAm = [66.2222 , 69 , 35.3333 , 38.2500];
% spikes_per_min_AATy = [293 , 498.5556 , 389.7778];

%for 3 max width
spikes_per_min_ctrl = [0 , 0 , 0.2222 ];
spikes_per_min_AAm = [0.1111 , 0 , 0.1111 , 5.5000];
spikes_per_min_AATy = [25 , 59.3333 , 33];

x = categorical({'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});
x = reordercats(x,{'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});
data = [mean(spikes_per_min_ctrl),mean(spikes_per_min_AAm),mean(spikes_per_min_AATy)]';
errhigh = [std(spikes_per_min_ctrl),std(spikes_per_min_AAm),std(spikes_per_min_AATy)];
errlow  = [std(spikes_per_min_ctrl),std(spikes_per_min_AAm),std(spikes_per_min_AATy)];

figure(254); subplot(1,3,2);
b = bar(x,data);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.14 0.14 0.8];
    b.CData(2,:) = [0.8 0.14 0.14];
    b.CData(3,:) = [0.14 0.8 0.14];
ylabel("Spikes Per Min");

hold on
    er = errorbar(x,data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
hold off

[h,p] = ttest2(spikes_per_min_ctrl , spikes_per_min_AAm , 'Vartype' , 'unequal')
[h,p] = ttest2(spikes_per_min_ctrl , spikes_per_min_AATy , 'Vartype' , 'unequal')
[h,p] = ttest2(spikes_per_min_AAm , spikes_per_min_AATy , 'Vartype' , 'unequal')
%% BAR GRAPH - product percentage of active cells * spikes per min

prod_ctrl = spikes_per_min_ctrl .* pct_spikes_ctrl; 
prod_AAm = spikes_per_min_AAm .* pct_spikes_AAm;
prod_AATy = spikes_per_min_AATy .* pct_spikes_AATy;

x = categorical({'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});
x = reordercats(x,{'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});
data = [mean(prod_ctrl),mean(prod_AAm),mean(prod_AATy)]';
errhigh = [std(prod_ctrl),std(prod_AAm),std(prod_AATy)];
errlow  = [std(prod_ctrl),std(prod_AAm),std(prod_AATy)];

figure(254); subplot(1,3,3);
b = bar(x,data);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.14 0.14 0.8];
    b.CData(2,:) = [0.8 0.14 0.14];
    b.CData(3,:) = [0.14 0.8 0.14];
ylabel("Product");

hold on
    er = errorbar(x,data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
hold off

[h,p] = ttest2(prod_ctrl , prod_AAm , 'Vartype' , 'unequal')
[h,p] = ttest2(prod_ctrl , prod_AATy , 'Vartype' , 'unequal')
[h,p] = ttest2(prod_AAm , prod_AATy , 'Vartype' , 'unequal')
%% BAR GRAPHS - pct of spikes in active cells 

%for 30s max width
% tot_spikes_ctrl = [10/4 , 200/9 , 109/4 ];
% tot_spikes_active_ctrl = [0/4 , 169/9 , 44/4]; 
% tot_spikes_AAm = [666/9 , 842/9 , 468/9 , 184/4];
% tot_spikes_active_AAm = [548/9 , 675/9 , 293/9 , 123/4];
% tot_spikes_AATy = [2864/9 , 4793/9 , 3825/9 ];
% tot_spikes_active_AATy = [2739/9 , 4693/9 , 3744/9];

%for 10s max width
% tot_spikes_ctrl = [0/4 , 49/4 , 143/9 ];
% tot_spikes_active_ctrl = [0/4 , 23/4 , 124/9]; 
% tot_spikes_AAm = [596/9 , 621/9 , 318/9 , 153/4];
% tot_spikes_active_AAm = [499/9 , 482/9 , 199/9 , 105/4];
% tot_spikes_AATy = [2637/9 , 4487/9 , 3508/9 ];
% tot_spikes_active_AATy = [2503/9 , 4381/9 , 3392/9];

%for 10s max width
tot_spikes_ctrl = [0/4 , 0/4 , 2/9 ];
tot_spikes_active_ctrl = [0/4 , 0/4 , 0/9]; 
tot_spikes_AAm = [1/9 , 0/9 , 1/9 , 22/4];
tot_spikes_active_AAm = [0/9 , 0/9 , 0/9 , 18/4];
tot_spikes_AATy = [225/9 , 534/9 , 297/9 ];
tot_spikes_active_AATy = [141/9 , 386/9 , 210/9];

x = categorical({'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});
x = reordercats(x,{'Control','Actin Arrest in Medium','Actin Arrest in TyNa'});

data = [tot_spikes_active_ctrl/tot_spikes_ctrl,tot_spikes_active_AAm/tot_spikes_AAm, tot_spikes_active_AATy/tot_spikes_AATy]
err = [std([tot_spikes_active_ctrl(1)/tot_spikes_ctrl(1),tot_spikes_active_ctrl(2)/tot_spikes_ctrl(2),tot_spikes_active_ctrl(3)/tot_spikes_ctrl(3)]) , ...  
    std([tot_spikes_active_AAm(1)/tot_spikes_AAm(1) , tot_spikes_active_AAm(2)/tot_spikes_AAm(2) , tot_spikes_active_AAm(3)/tot_spikes_AAm(3) , tot_spikes_active_AAm(4)/tot_spikes_AAm(4)]) , ...
    std([tot_spikes_active_AATy(1)/tot_spikes_AATy(1) , tot_spikes_active_AATy(2)/tot_spikes_AATy(2) , tot_spikes_active_AATy(3)/tot_spikes_AATy(3)]) ];

figure(265482);
b=bar(x,data); 
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.14 0.14 0.8];
    b.CData(2,:) = [0.8 0.14 0.14];
    b.CData(3,:) = [0.14 0.8 0.14];
ylabel("Percentage of Spikes from Active Cells");

hold on
    er = errorbar(x,data,err,err);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
hold off

%% Some sort of frequency analysis??

Fs = fps;               % Sampling frequency in Hz                 
T = 1/Fs;             % Sampling period       
L = size(DFF,2);             % Length of signal in milliseconds 
t = (0:L-1)*T;        % Time vector

%%
y = fft(DFF(196,:));

P2 = abs(y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure(3452);
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%%

n = length(zen(270,:));
fs = fps;
t = (0:n-1)/fs;

for n = 1:size(DFF,1)
y = fft(DFF(n,:));
power = (abs(y).^2)/n;
FFTkymo(n,:) = power;
end

y = fft(DFF(177,:)); 

power = (abs(y).^2)/n;
f = (0:n-1)*(fs/n);
figure(2352);
% plot(f,power);
plot(f(1:floor(n/2)),power(1:floor(n/2)))
hold on;

%%
figure(25);
imagesc(FFTkymo(1:floor(n/2)));




