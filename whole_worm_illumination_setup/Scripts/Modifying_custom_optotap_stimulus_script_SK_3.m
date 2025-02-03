%%%% This is the final version of the code updated on 07/25/2020
%%%% variable light_before_tap is applicable when we are giving variable
%%%% amount of light stim before tap. LED_before_tap is more universal

close all
clear
clc

% plot behavior profiles under different LED intensities
cd('/projects/LEIFER/Sandeep/Data/20221115_array_light_mix_tap_duration_details_Sandeep_Adtran_KP4_glr-1_noATR_day2')
% cd ('/projects/LEIFER/Sandeep/Data/20220328_array_light_mix_tap_duration_details_Sandeep_CGZ195/Data20220328_175514')
% cd ('D:\Expt_Data\AML358_ACR2_data\With_tap\AML358_20200306') %%% the file ending with 134311 has 0.5 only light and 0.5 only Tap
% cd('/projects/LEIFER/Sandeep/Data/20200306_array_light_mix_tap_duration_details_AML358_Sandeep')

%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%
analysis_workspace_data=11;
analysis_new_data=1;
current_rails_dur=14*3;  %%% 14 fps x time in sec.

if analysis_workspace_data==1
   [file,path] = uigetfile('*.mat');
   fullpath=fullfile(path,file);
   disp(fullpath);
   load(fullpath)
   analysis_new_data=0;
%    load('/projects/LEIFER/Sandeep/Data/20200109_array_light_mix_tap_details_AKS337.9.6.b_Sandeep/20200109_workspace_data.mat') 
end

plot_all_beh_ratios=1;
plot_various_int=1;
plot_comparison_reverse_beh=11;

if analysis_new_data==1
% plot behavior profiles under different LED intensities

load('reference_embedding_3behaviors.mat')
relevant_track_fields = {'BehavioralTransition','Frames','TapVoltages','LEDVoltages','Velocity'};

%%%%%%% select folders
folders_optotap = getfoldersGUI();

%%%%% For waitbar
% % % % waitbar_var = 0;
% % % % disp_waitbar = waitbar(0, 'Please wait...');

% specify tracks of stimuli condition
LED_and_Tap=0;
all_Tap=0; % includes zeros LED power
LED_only=1;
LED_only_varying_stim=11;
Tap_only=1;
Tap_varying_opto_stim=0;
LED_before_tap=0; %%%% e.g. 3, 7, or 14 (when using Tap_varying_opto_stim)
read_from_text_file=1; %%% use it as 1 when you want to read some info directly from parameters.txt file (e.g. 03/06/2020 data)

%Global calibration uses parameters.csv for power conversion
%otherwise uses the "parameters.txt" in each folder
global_calibration=0;

%load stimuli.txt from the first experiment
num_stimuli = 1;
time_window_before = 140;
time_window_after = 140;
total_window_frames = time_window_before+time_window_after+1;
fps = 14;
stim_similarity_thresh = 1.01;  % 1.11 or 1.2

number_of_behaviors = max(L(:)-1);
stimulus_intensities = [];
all_behavior_transitions_for_frame = {};
all_behavior_annotations_for_frame = {};

%% behavioral rate compare
for folder_index = 1:length(folders_optotap)
    sprintf('Analyzing data %d out of %d datasets', folder_index, length(folders_optotap))

% % % %     if ~ishghandle(disp_waitbar)
% % % %     error('User closed waitbar.');
% % % %     end
% % % %     waitbar_var = waitbar_var + 1;
% % % %     waitbar(waitbar_var / length(folders_optotap), disp_waitbar, sprintf('%d of %d', waitbar_var, length(folders_optotap)));
% % % %        
    %load the tracks for this folder
    [folder_tracks, folder_indecies_revstim_ret, track_indecies_revstim_ret] = loadtracks(folders_optotap{folder_index},relevant_track_fields);
    folder_tracks = BehavioralTransitionToBehavioralAnnotation(folder_tracks);
    %generate the Behavior matricies
    folder_tracks = get_behavior_triggers(folder_tracks);
    
    % separate the tracks based on stimuli conditions
    % [sti_tracks,TAPonly_tracks,LEDonly_tracks]=get_tracks_with_sti(folder_tracks);
    %specify which conditions to analyze
    current_tracks=folder_tracks;
    
    current_param = load_parameters(folders_optotap{folder_index});
    LEDVoltages = load([folders_optotap{folder_index}, filesep, 'LEDVoltages.txt']);
    LEDVoltages=abs(LEDVoltages);  %%% I am doing this in the case when the control intensty is zero (by giving negative stim)

    TapVoltages = load([folders_optotap{folder_index}, filesep, 'TapVoltages.txt']);
    tap_parameters= importdata([folders_optotap{folder_index}, filesep, 'parameters.txt']);
    
    %convert LEDVoltages to power using global calibration
    if global_calibration
        LEDPowers = round(LEDVoltages .* current_param.avgPower500 ./ 5);
    else % Use the parameters.txt (LabView input data) in each folder
        folder_paras=readtable([folders_optotap{folder_index}, filesep, 'parameters.txt'],'Delimiter','\t');
% %         LEDPowers = round(LEDVoltages.*folder_paras.VoltageToPower); %%% Mochi's code
        LEDPowers = ceil(LEDVoltages.*folder_paras.VoltageToPower); %% I am not rounding or else the control intensity gets mixed with zero
        LEDPowers(end)=0; % correct for the last peak, otherwise undetected
        %find when each stimuli is played back by convolving the time
        %reversed stimulus (cross-correlation)
    end
    
    LED_switching_points=[0 diff(LEDPowers)];

% % % %     [~, LED_peaks_locs,LED_peaks_widths,~] =
% findpeaks(LED_switching_points, 'MinPeakDistance',14); %%% original code
    [LED_peaks_amplitude, LED_peaks_locs,LED_peaks_widths,~] = findpeaks(LEDPowers, 'MinPeakDistance',14);
    [~, Tap_peaks] = findpeaks(TapVoltages, 'MinPeakDistance',14);
    
    max_allowed_stim_duration_diff=5;                 
    raw_stim_peaks_locs=[];
    raw_stim_peaks_height=[];
    raw_stim_duration=[];
    for i=1:size(LED_peaks_widths,2)
        if (abs(current_rails_dur - LED_peaks_widths(1,i)))<max_allowed_stim_duration_diff
        raw_stim_peaks_locs=[raw_stim_peaks_locs LED_peaks_locs(1,i)];  %%% these are all the peaks with same width
        raw_stim_peaks_height=[raw_stim_peaks_height LED_peaks_amplitude(1,i)];
        raw_stim_duration=[raw_stim_duration LED_peaks_widths(1,i)];
        end
    end

    stim_of_certain_width=zeros(size(LEDPowers))';
    for i=1:size(raw_stim_peaks_locs,2)
        stim_of_certain_width([raw_stim_peaks_locs(1,i):(raw_stim_peaks_locs(1,i)+raw_stim_duration(1,i))],1)=raw_stim_peaks_height(1,i);
    end
        
    %%%%%% find stimuli time point based on conditions
    if Tap_varying_opto_stim==1
    for i=1:size(LED_peaks_locs,2)
       time_diff(i)=min(abs(Tap_peaks-LED_peaks_locs(i))); 
    end

    for i=1:size(light_before_tap,2)
        cluster_assignments{:,i}= find(light_before_tap(:,i)==time_diff);
    end

    for m=1:size(light_before_tap,2)
        peak_locations_time_diff= LED_peaks_locs(cluster_assignments{:,m});    %%%%% use this for LED_and_tap
        all_light_tap_peaks_array{:,m}=peak_locations_time_diff+light_before_tap(:,m);  %%% throught this we will find all_tap
        all_light_tap_peaks_cluster{:,m}=all_light_tap_peaks_array{:,m};
    end

    all_light_tap_peaks_array=cell2mat(all_light_tap_peaks_array);
    tap_only_peaks=setdiff(Tap_peaks,all_light_tap_peaks_array);

    ab=find(light_before_tap==LED_before_tap);
    end

    if LED_and_Tap==1
        peak_locations=all_light_tap_peaks_cluster{:,ab};  %%%  here ab is the index of the light duraton before tap   
    elseif all_Tap==1
        peak_locations=Tap_peaks;
    elseif LED_only==1
        peak_locations=LED_peaks_locs;
    elseif LED_only_varying_stim==1
        peak_locations=raw_stim_peaks_locs;
    elseif Tap_only==1
        peak_locations=tap_only_peaks;
    elseif Tap_varying_opto_stim==1
        peak_locations=union(all_light_tap_peaks_cluster{:,ab},tap_only_peaks);
        correspond_LED_power=LEDPowers(peak_locations);
    end
    
%     %%
%     figure;
%     x_axis=linspace(1/14,1800/60,25200);
%     plot(LEDPowers,'-r','LineWidth',1.2)
%     hold on
%     plot(TapVoltages,'-b','LineWidth',1.2)
%     hold on;
%     plot(peak_locations,LEDPowers(peak_locations),'og')
%     xlabel('Time (sec)','fontsize',14); ylabel('Stimulus intensity (\muW/mm^2)','fontsize',14);   %%%% Writing x and y label. The fontsize of the label is also mentioned.                                                         %%%% Writing the title of the plot
% %     xlim([0 250]); 
%     ylim([0 55]) 
%     set(gca,'fontsize',14)    
%     box off
    %%
    %loop through the peaks to update light powers
    for peak_index = 1:length(peak_locations)
%         disp(peak_index)
        
        if Tap_only
%             sprintf('Tap_only')
            stimulus_intensities=0;
            current_stim_index=1;
            all_behavior_transitions_for_frame{1} = cell(1,total_window_frames);
            all_behavior_annotations_for_frame{1} = cell(1,total_window_frames);
            
        else  %get the stimulus intensity for this peak
%             sprintf('not Tap_only')
            current_stim_power = LEDPowers(peak_locations(peak_index));
            if current_stim_power==0
%                 sprintf('current_stim_power==0')
                if ismember(0,stimulus_intensities) % 0 already exists
                    current_stim_index=find(stimulus_intensities==0);
                else
                    stimulus_intensities = [stimulus_intensities,current_stim_power];
                    current_stim_index = length(stimulus_intensities);
                    all_behavior_transitions_for_frame{current_stim_index} = cell(1,total_window_frames);
                    all_behavior_annotations_for_frame{current_stim_index} = cell(1,total_window_frames);
                end
            else %non zero LED power
%                 sprintf('current_stim_power=~0')
                %see if this stim_power already exists
                current_stim_index = find(and(stimulus_intensities < current_stim_power*stim_similarity_thresh,stimulus_intensities > current_stim_power/stim_similarity_thresh));
%                 disp(stimulus_intensities)
%                 disp(current_stim_power*stim_similarity_thresh)
%                 disp(current_stim_power/stim_similarity_thresh)
                
                if isempty(current_stim_index)
%                     sprintf('current_stim_index_is)empty')
                    %no entry yet
                    stimulus_intensities = [stimulus_intensities,current_stim_power];
                    current_stim_index = length(stimulus_intensities);
                    all_behavior_transitions_for_frame{current_stim_index} = cell(1,total_window_frames);
                    all_behavior_annotations_for_frame{current_stim_index} = cell(1,total_window_frames);
                end
            end
        end
        
        %for every time a stimulus is delivered, look at a certain range of frames
        for frame_shift = -time_window_before:time_window_after
            current_frame = peak_locations(peak_index) +frame_shift;
            if current_frame <= length(LEDPowers) && current_frame >= 1
                %make sure the current frame is in range
                %cut up tracks to each frame
                tracks_on_critical_frame = FilterTracksByTime(current_tracks,current_frame, current_frame);
                if ~isempty(tracks_on_critical_frame)
                    all_behavior_transitions_for_frame{current_stim_index}{frame_shift+time_window_before+1} = [all_behavior_transitions_for_frame{current_stim_index}{frame_shift+time_window_before+1}, tracks_on_critical_frame.Behaviors];
                    all_behavior_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1} = [all_behavior_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1}, tracks_on_critical_frame.BehavioralAnnotation];
                end
            end
        end
    end
    clear all_light_tap_peaks_array all_light_tap_peaks_cluster
end

%sort the stimulus intensities
[stimulus_intensities, sort_index] = sort(stimulus_intensities);
all_behavior_transitions_for_frame = all_behavior_transitions_for_frame(sort_index);
all_behavior_annotations_for_frame= all_behavior_annotations_for_frame(sort_index);

% %% 1 plot the transition rates as a function of time
% for stimulus_index = 1:length(stimulus_intensities)
%     % plot the transition rates centered on stim delivery
%     transition_counts_for_frame = zeros(number_of_behaviors,total_window_frames);
%     transition_rate_for_frame = zeros(number_of_behaviors,total_window_frames);
%     transition_std_for_frame = zeros(number_of_behaviors,total_window_frames);
%     for frame_index = 1:total_window_frames
%         transitions_for_frame = all_behavior_transitions_for_frame{stimulus_index}{frame_index};
%         transition_counts_for_frame(:,frame_index) = sum(transitions_for_frame,2);
%         transition_rate_for_frame(:,frame_index) = sum(transitions_for_frame,2)./size(transitions_for_frame,2).*fps.*60;
%         transition_std_for_frame(:,frame_index) = sqrt(sum(transitions_for_frame,2))./size(transitions_for_frame,2).*fps.*60;
%     end
%
%     track_n = round(mean(arrayfun(@(x) size(x{1},2), [all_behavior_transitions_for_frame{stimulus_index}])));
%     my_colors = behavior_colors;
%
%     figure
%     hold on
%     for behavior_index = 1:number_of_behaviors
%         transition_n = sum(transition_counts_for_frame(behavior_index,:));
%         plot(-time_window_before/fps:1/fps:time_window_after/fps, transition_rate_for_frame(behavior_index,:), '-', 'color', my_colors(behavior_index,:),'Linewidth', 3,'DisplayName',[behavior_names{behavior_index}, ' (', num2str(transition_n),' transitions)']);
%     end
%     hold off
%     xlabel('Time (s)') % x-axis label
%     ylabel('Transition Rate (transitions/min)') % y-axis label
%     title(['Stimulus Intensity = ', num2str(stimulus_intensities(stimulus_index)), ' (n = ', num2str(track_n), ' tracks)']);
%     legend('show');
%     ax = gca;
%     ax.FontSize = 10;
%     axis([-10 10 0 35])
% end
n_sti=length(stimulus_intensities);
behavior_counts_for_frame = zeros(number_of_behaviors,n_sti,total_window_frames);
behavior_ratios_for_frame = zeros(number_of_behaviors,n_sti,total_window_frames);

for stimulus_index = 1:n_sti
    % plot the transition rates centered on stim delivery
    total_counts_for_frame = zeros(1,total_window_frames);
    for frame_index = 1:total_window_frames
        for behavior_index = 1:number_of_behaviors
            behavior_counts_for_frame(behavior_index,stimulus_index,frame_index) = sum(all_behavior_annotations_for_frame{stimulus_index}{frame_index}==behavior_index);
        end
        behavior_ratios_for_frame(:,stimulus_index,frame_index) = behavior_counts_for_frame(:,stimulus_index,frame_index)./sum(behavior_counts_for_frame(:,stimulus_index,frame_index)); %get ratio
    end
end
end

%%% 2 plot the behavioral ratios as a function of time
%%
if plot_all_beh_ratios==1
n_tracks=zeros(1,n_sti); %number of tracks in each sti intensities
for stimulus_index = 1:n_sti
    
    n_tracks(stimulus_index) = round(mean(arrayfun(@(x) size(x{1},2), [all_behavior_transitions_for_frame{1,stimulus_index}])));
    my_colors = behavior_colors;
    figure
    light_pulse = area([0 current_rails_dur/fps], [1 1],'facecolor',[1 0 0],'LineStyle','none');
    hold on;grid off
    for behavior_index = 1:number_of_behaviors
        plot(-time_window_before/fps:1/fps:time_window_after/fps, squeeze(behavior_ratios_for_frame(behavior_index,stimulus_index,:)), '-', 'color', my_colors(behavior_index,:),'Linewidth', 3,'DisplayName',behavior_names{behavior_index});
    end
       hold off
    xlabel('Time (s)') % x-axis label
    ylabel('Behavioral Ratio') % y-axis label
%     title(['Strain: F; Stimulus Intensity = ', num2str(stimulus_intensities(stimulus_index)),'\muW/mm^2' ' (n = ', num2str(n_tracks(stimulus_index)), ' tracks)']);
    title(['CGZ221', ' (n = ', num2str(n_tracks(stimulus_index)), ' tracks)']);
%     tap_line=xline(0,'k--');
   % set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(tap_line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('show','Location','northwest');
    ax = gca;
    ax.FontSize = 13;
    axis([-time_window_before/fps time_window_after/fps 0 1])
    alpha(0.1)
    box off;
%     xlim([-7 3]); 
%     ylim([0 0.3])
end
end
%% plot the behavioral ratios for various intensities on the same plot
if plot_various_int==1
for behavior_index = [2,3] %forward3; forward 5; fast reverse; turns
    my_colors = lines(n_sti);
    figure;
    light_pulse = area([0 current_rails_dur/fps], [1 1],'facecolor',[1 0 0],'LineStyle','none');
    hold on; grid on
    for stimulus_index = 1:n_sti
        %track_n = round(mean(arrayfun(@(x) size(x{1},2), [all_behavior_transitions_for_frame{1,stimulus_index}])));
        plot(-time_window_before/fps:1/fps:time_window_after/fps, squeeze(behavior_ratios_for_frame(behavior_index,stimulus_index,:)), '-', 'color', my_colors(stimulus_index,:),'Linewidth', 3,...
            'DisplayName',[num2str(stimulus_intensities(stimulus_index)), 'uW/mm2 (n = ', num2str(n_tracks(stimulus_index)),')']);
    end
    hold off
    xlabel('Time (s)') % x-axis label
    ylabel('Behavioral Ratio') % y-axis label
    title(['Strain name: CGZ221 ', 'Beh: ' behavior_names{behavior_index}]);
    tap_line=line(0,'k--');
    set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(tap_line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ax = gca;
    ax.FontSize = 14;
    legend('show','Location','northwest');
    ylim([0 0.5])
    alpha(0.1)
    box off; grid off;
end
end
%% plot how we picked the GWN range

if plot_comparison_reverse_beh==1
%optotap_behavioral_ratio_percent_changes = percent_change_above_baseline(squeeze(behavior_ratios_for_frame(:,4,:)));
behavior_index = 3;
% [behavior_percent_change, behavior_baselines, behavior_max, behavior_min] = percent_change_above_baseline(squeeze(behavior_ratios_for_frame(behavior_index,:,:)));
% y = behavior_max';
x = stimulus_intensities;
err_boot=zeros(1,length(x));
% standard error of the mean
sem=zeros(1,length(x));
% find the max ratio within the time window
% % % % [y2,mi]= max(behavior_ratios_for_frame(mn,:,301:end),[],3);  %%%original code
[y2,mi]= max(behavior_ratios_for_frame(behavior_index,:,:),[],3);

for stimulus_index = 1:n_sti
    % create binary data at the frame when the ratio of behavior is max
% % %
% binary_data=(all_behavior_annotations_for_frame{stimulus_index}{mi(stimulus_index)}==behavior_index);%% original code
    binary_data=(all_behavior_annotations_for_frame{stimulus_index}{mi(stimulus_index)}==behavior_index);  %%%% Adding the same index number is very important
        
    % get standard error from bootstraped data
    % bootratio=bootstrp(20000,@mean,binary_data);
    % err_boot(stimulus_index)=std(bootratio);
    % mathematical formula of standard error of the mean
    sem(stimulus_index)=std(binary_data)/sqrt(length(binary_data));
end
% figure('Position',[100,100,800,600])
figure;
errorbar(x, y2,sem, 'bo-', 'LineWidth',2,'Markersize',10)
for stimulus_index = 1:n_sti
    %track_n = round(mean(arrayfun(@(x) size(x{1},2), [all_behavior_transitions_for_frame{1,stimulus_index}])));
    text(x(stimulus_index), y2(stimulus_index), ['   n=', num2str(n_tracks(stimulus_index))]);
end
ax = gca;
ax.XTick = stimulus_intensities;
% ax.YTick = [0 0.3 0.6];
ax.FontSize = 13;

xlabel('Stimulus Intensity (\muW/mm^2)') % x-axis label
ylabel('Fast Reverse Ratio') % y-axis label
title('10 \mul ATR in the tap assay setup');
ylim([0,0.5]);
% xlim([0 60])
% axis tight;
box off
end
%%%%%
% F=findall(0,'type','figure','tag','TMWWaitbar');
% delete(F)

%%%
% % % % figure;
% % % % x_axis=linspace(1/14,1800/60,25200);
% % % % plot(x_axis,LEDPowers,'-r','LineWidth',1.2)
% % % % hold on
% % % % plot(x_axis,TapVoltages*10,'-b','LineWidth',1.2)
% % % % xlabel('Time (sec)','fontsize',14); ylabel('Stimulus intensity (\muW/mm^2)','fontsize',14);   %%%% Writing x and y label. The fontsize of the label is also mentioned.                                                         %%%% Writing the title of the plot
% % % % % xlim([0 250]); ylim([0 55]) 
% % % % set(gca,'fontsize',14)    
% % % % box off

%%% Clear some variables
% clear plot_all_beh_ratios plot_various_int plot_comparison_reverse_beh analysis_workspace_data analysis_new_data