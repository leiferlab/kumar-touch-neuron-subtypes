%%%% This code turns LED ON for different time durations before the tap
clear
clc
close all

% powers are for LED light; voltages are for tapping
% sti is the length of stimuli (ligh is on)
% The first power input is regarded as default for non-tap 
dur=25200;
period=840;
prop_no_tap=0;
prop_no_light=0.2;
power_array=[50];
light_stim_duration_array=[42 35 31];
stim_duration_after_tap=28;               %%%% standard value of 112
taps=1;
%%%%%%%% hold_before_tap variable is deleted in this version of the code
powers=zeros(1,dur);
voltages=zeros(1,dur);
vol_tap=5;
cycle=floor(dur/period);
number_tap=round((1-prop_no_tap)*cycle);
number_light=round((1-prop_no_light)*cycle);
number_both=round((1-prop_no_light-prop_no_tap)*cycle);
number_only_light=number_light-number_both;

%%% modify the power_array size to match number_both
p_len=length(power_array);
n=floor(number_both/p_len);
m=power_array(randperm(p_len,mod(number_both,p_len)));
new_power_array=[repmat(power_array,1,n),m];

%%%%%%%%%% modify the hold before tap array size 
hold_tap_len=length(light_stim_duration_array);
new_n=floor(number_both/hold_tap_len);
new_m=light_stim_duration_array(randperm(hold_tap_len,mod(number_both,hold_tap_len)));
new_hold_tap_array=[repmat(light_stim_duration_array,1,new_n),new_m];
random_hold_tap_array=new_hold_tap_array;

for_light_only_array=repelem(light_stim_duration_array,ceil(cycle/length(light_stim_duration_array)));
random_hold_tap_array_light_only=for_light_only_array(randperm(cycle));

%%%%%%% power and voltage in one period
tap_int=(stim_duration_after_tap)/(taps); %tap intervals
p_per=zeros(1,period);
v_per=zeros(1,period);

ini_sti=period-stim_duration_after_tap+1;
for i =1:taps % tap_voltages in one standard period
    v_per(ini_sti+round(tap_int*(i-1)))=vol_tap;
end

%%%%%%%% randomly select which cycles to apply conditions
rand_ind_cycle=randperm(cycle);
light_only_index=rand_ind_cycle(1:number_only_light);
both_index=rand_ind_cycle(number_only_light+1:number_light);
tap_index=rand_ind_cycle(cycle-number_tap+1:end);

for i=1:number_tap
    cycle_inx_tap=tap_index(i);
    end_frame=cycle_inx_tap*period;
    ini_frame=end_frame-period+1;
    voltages(ini_frame:end_frame)=v_per;
end

for i=1:number_both
    
    cycle_inx_both=tap_index(i);
    end_frame=cycle_inx_both*period;
    ini_frame=end_frame-period+1;
    voltages(ini_frame:end_frame)=v_per;
        
    ini_sti=period-random_hold_tap_array(i)+1;
    p_per(ini_sti:period)=1;
        
    cycle_inx_both=both_index(i);
    power=new_power_array(i);
    end_frame=cycle_inx_both*period;
    ini_frame=end_frame-period+1;
    powers(ini_frame:end_frame)=p_per.*power;
    p_per=zeros(1,period);
end

%%%%%% for light only cycles, we use the first power input

%%%%%% for_light_only=[repmat(power_array,1, round(cycle/p_len)),[]];
for_light_only=repelem(power_array,ceil(cycle/p_len));
power_array_light_only=for_light_only(randperm(cycle));

for i=1:number_only_light
    cycle_inx_light=light_only_index(i);
    power=power_array_light_only(i);
    end_frame=cycle_inx_light*period;
    ini_frame=end_frame-period+1;
    
    ini_sti=period-random_hold_tap_array_light_only(i)+1;
    p_per(ini_sti:period)=1;
    powers(ini_frame:end_frame)=p_per.*power;
    p_per=zeros(1,period);
end

% figure; plot(linspace(1/14,1800,25200),powers,'-r');
% hold on;
% plot(linspace(1/14,1800,25200),voltages*10,'-b');
% % xlim([0 300])

%%
%%%%%%%% To detect cycle parameters
all_light_duration=light_stim_duration_array;
light_before_tap=all_light_duration-stim_duration_after_tap;

LEDPowers = powers;
TapVoltages= voltages;

LED_and_Tap=0;
all_Tap=0; % includes zeros LED power
LED_only=0;
Tap_only=0;
Tap_varying_opto_stim=1;
LED_before_tap=7; 
normalized_stimuli = 1;

LED_switching_points=[0 diff(LEDPowers)];

[~, LED_peaks] = findpeaks(LED_switching_points, 'MinPeakDistance',14);
[~, Tap_peaks] = findpeaks(TapVoltages, 'MinPeakDistance',14);

%%%%%% find stimuli time point based on conditions
for i=1:size(LED_peaks,2)
   time_diff(i)=min(abs(Tap_peaks-LED_peaks(i))); 
end

for i=1:size(light_before_tap,2)
    cluster_assignments{:,i}= find(light_before_tap(:,i)==time_diff);
end

for m=1:size(light_before_tap,2)
    peak_locations_time_diff{:,m}= LED_peaks(cluster_assignments{:,m});    %%%%% use this for LED_and_tap
    all_light_tap_peaks_array{:,m}=peak_locations_time_diff{:,m}+light_before_tap(:,m);  %%% throught this we will find all_tap
    all_light_tap_peaks_cluster{:,m}=all_light_tap_peaks_array{:,m};
end

all_light_tap_peaks_array=cell2mat(all_light_tap_peaks_array);
tap_only_peaks=setdiff(Tap_peaks,all_light_tap_peaks_array);

ab=find(light_before_tap==LED_before_tap);

if LED_and_Tap
    peak_locations=all_light_tap_peaks_cluster{:,ab}    %%%  here ab is the index of the light duraton before tap   
elseif all_Tap
    peak_locations=Tap_peaks
elseif LED_only
    peak_locations=LED_peaks
elseif Tap_only
    peak_locations=tap_only_peaks
elseif Tap_varying_opto_stim
    peak_locations=union(all_light_tap_peaks_cluster{:,ab},tap_only_peaks)    
end

%%
figure;
x_axis=linspace(1/14,1800/60,25200);
plot(LEDPowers,'-r','LineWidth',1.2)
hold on
plot(TapVoltages,'-b','LineWidth',1.2)
hold on;
plot(peak_locations,LEDPowers(peak_locations),'og')
xlabel('Time (sec)','fontsize',14); ylabel('Stimulus intensity (\muW/mm^2)','fontsize',14);   %%%% Writing x and y label. The fontsize of the label is also mentioned.                                                         %%%% Writing the title of the plot
%     xlim([0 250]); 
ylim([0 155]) 
set(gca,'fontsize',14)    
box off