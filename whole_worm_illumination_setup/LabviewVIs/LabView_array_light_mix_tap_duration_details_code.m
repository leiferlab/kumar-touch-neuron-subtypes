%%%% This is the new script
%light_stim_duration_array=[42 70 98];
stim_duration_after_tap=28;
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