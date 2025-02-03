%%%% This code will analyze just opto stim experiments in the tap assay setup

close all
clear
clc

addpath(genpath('/projects/LEIFER/Sandeep/github/leifer-Behavior-Triggered-Averaging-Tracker'));

cd ('/projects/LEIFER/Sandeep/Data/20220407_array_light_mix_tap_duration_details_Sandeep_CGZ546');

%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%
analysis_workspace_data=11;
analysis_new_data=1;
plot_each_folder=11;
plot_all_intensities_together=1;
plot_all_behavior_together=1;
plot_behavior_heatmaps=1;
plot_barplots=1;
plot_velocity_traces=1;
plot_velocity_heatmaps=1;
stim_duration=3;
current_rails_dur=14*stim_duration;  %%% 14 fps x time in sec.
ellipse_ratio_threshold=0.8;
time_window_before = 140;
time_window_after = 140;
total_window_frames = time_window_before+time_window_after+1;
fps = 14;
stim_intensity_max=250;
stim_intensity_min=150;

if analysis_workspace_data==1
   [file,path] = uigetfile('*.mat');
   fullpath=fullfile(path,file);
   disp(fullpath);
   load(fullpath)
   analysis_new_data=0;
%    load('/projects/LEIFER/Sandeep/Data/20200109_array_light_mix_tap_details_AKS337.9.6.b_Sandeep/20200109_workspace_data.mat') 
end

if analysis_new_data==1
% plot behavior profiles under different LED intensities

load('reference_embedding_3behaviors.mat')
relevant_track_fields = {'BehavioralTransition','Frames','LEDVoltages','Velocity','Eccentricity','Centerlines','Path'};

%%%%%%% select folders
folders_optotap = getfoldersGUI();
stimulus_intensities = [];

%% behavioral rate compare
behavior_state_stim_folder=cell(1,length(folders_optotap));
velocity_state_stim_folder=cell(1,length(folders_optotap));

for folder_index = 1:length(folders_optotap)
    sprintf('Analyzing data %d out of %d datasets', folder_index, length(folders_optotap))
 
    %load the tracks for this folder
    [folder_tracks, ~, ~] = loadtracks(folders_optotap{folder_index},relevant_track_fields);
    folder_tracks = BehavioralTransitionToBehavioralAnnotation(folder_tracks);
    %generate the Behavior matricies
    current_tracks = get_behavior_triggers(folder_tracks);
    
    current_param = load_parameters(folders_optotap{folder_index});
    LEDVoltages = load([folders_optotap{folder_index}, filesep, 'LEDVoltages.txt']);
    LEDVoltages=abs(LEDVoltages);  %%% I am doing this in the case when the control intensty is zero (by giving negative stim)

    % convert LEDVoltages to power 
    
    folder_paras=readtable([folders_optotap{folder_index}, filesep, 'parameters.txt'],'Delimiter','\t');
    LEDPowers = ceil(LEDVoltages.*folder_paras.VoltageToPower); %% I am not rounding or else the control intensity gets mixed with zero
    LEDPowers(end)=0; % correct for the last peak, otherwise undetected

    [LED_peaks_amplitude, LED_peaks_locs,LED_peaks_widths,~] = findpeaks(LEDPowers, 'MinPeakDistance',14);
    
    stim_intensities=unique(LED_peaks_amplitude);
    behavior_state_stim=cell(size(stim_intensities,2),1);
    velocity_state_stim=cell(size(stim_intensities,2),1);
    
    size_current_track_matrix(folder_index,1)=size(current_tracks,2);
    
    %%%% if you want to get rid of folders with very few tracks
% % % %     if size(current_tracks,2)<=50
% % % %         behavior_state_stim_folder{:,folder_index}={};
% % % %         continue
% % % %     end
    
    for stim_inten=1:size(stim_intensities,2)
        behavior_states_in_track=[];
        velocity_states_in_track=[];
        locs_track_index_array=[];
        
        for track_index=1:size(current_tracks,2)
            
            voltage_track_index=current_tracks(track_index).LEDVoltages;
            LEDPowers_track_index = ceil(voltage_track_index.*folder_paras.VoltageToPower);
            
            %%% find peaks on LED power of each track
            [LEDPowers_track_index_amplitude, LEDPowers_track_index_locs,LEDPowers_track_index_widths,~] = findpeaks(LEDPowers_track_index, 'MinPeakDistance',14);
            stim_index=(LEDPowers_track_index_amplitude==stim_intensities(stim_inten));
            locs_index=LEDPowers_track_index_locs(stim_index);
            
            if isempty(locs_index)
                continue
            end
            
            %%%% if pipeline does'nt assign behavior to one of the three states then we will take velocity into consideration to assign forward or reverse
            original_ellipse_ratio_for_correction=current_tracks(track_index).Eccentricity;
            smooth_ellipse_ratio_for_correction=smooth(original_ellipse_ratio_for_correction,10)';
            
            original_velocity_for_correction=current_tracks(track_index).Velocity;
            smooth_velocity_for_correction=smooth(original_velocity_for_correction,10);   %%%% smoothing the velocity data
            current_tracks(track_index).Velocity=smooth_velocity_for_correction;  %%% redefinig the velocity info with smootherd data
            smooth_velocity_for_correction=smooth_velocity_for_correction';
            
            for xyz=1:size(current_tracks(track_index).Eccentricity,2)
                if smooth_ellipse_ratio_for_correction(1,xyz)<ellipse_ratio_threshold
                    current_tracks(track_index).BehavioralAnnotation(1,xyz)=2;
                end
            end

            original_data_pre_correction=current_tracks(track_index).BehavioralAnnotation;
            data_unassigned_pre_correction=find(original_data_pre_correction==0); %%%% when beh is unassigned

            corrected_data_for_missing_behaviors=original_data_pre_correction;
            for xx=[data_unassigned_pre_correction]
                if smooth_velocity_for_correction(1,xx)<=-0.05
                    corrected_data_for_missing_behaviors(1,xx)=3;
                else
                    corrected_data_for_missing_behaviors(1,xx)=1;   
                end
            end

            current_tracks(track_index).BehavioralAnnotation=corrected_data_for_missing_behaviors;  %%%% we are reassigning behavior index for missing frames

            for mn=1:size(locs_index,2)
               behavior_states_in_dummy=[];
               velocity_states_in_dummy=[];
                
               if locs_index(mn)<=time_window_before
                   continue
               end
               
               %%% if track ends inbetween stim then ignore those tracks
               if size(current_tracks(track_index).Frames,2)-locs_index(mn)<=stim_duration*fps  
                   continue
               end
               
               %%% here I am adding NaNs in case some tracks are not more than 10 sec longer since stim initiation 
               behavior_on_track=current_tracks(track_index).BehavioralAnnotation;
               
               if size(behavior_on_track,2)-locs_index(mn)<=time_window_after 
                   original_size=size(current_tracks(track_index).BehavioralAnnotation,2);
                   nans_to_add=time_window_after-(original_size-locs_index(mn));
                   behavior_on_track(1,original_size+1:original_size+nans_to_add)=missing;
                   smooth_velocity_for_correction(1,original_size+1:original_size+nans_to_add)=missing;
               end
               
%                %%% ignore the cases when stim occured during turns
%                if any(current_tracks(track_index).BehavioralAnnotation(1,(locs_index(mn)-1:locs_index(mn)+1))==2) 
%                    continue
%                end
%                
%                %%% ignore the cases when stim occured during reversals
%                if any(current_tracks(track_index).BehavioralAnnotation(1,(locs_index(mn)-1:locs_index(mn)+1))==3) 
%                    continue
%                end
               
               %%% arranging the data so that I can use it for plottng individual worm videos
               if stim_intensity_min<=LEDPowers_track_index(locs_index(mn)) && stim_intensity_max>=LEDPowers_track_index(locs_index(mn))
                   dummy_locs_track_index(:,1)=track_index;
                   dummy_locs_track_index(:,2)=locs_index(mn);
                   locs_track_index_array=[locs_track_index_array; dummy_locs_track_index];
               end
               
               behavior_states_in_dummy=[behavior_states_in_dummy; behavior_on_track(1,(locs_index(mn)-time_window_before):(locs_index(mn)+time_window_after))];
               velocity_states_in_dummy=[velocity_states_in_dummy; smooth_velocity_for_correction(1,(locs_index(mn)-time_window_before):(locs_index(mn)+time_window_after))];
  
            end
            
            behavior_states_in_track=[behavior_states_in_track; behavior_states_in_dummy];
            velocity_states_in_track=[velocity_states_in_track; velocity_states_in_dummy];
        end
        behavior_state_stim{stim_inten,:}=behavior_states_in_track;
        velocity_state_stim{stim_inten,:}=velocity_states_in_track;
    end
    behavior_state_stim_folder{:,folder_index}=behavior_state_stim;
    velocity_state_stim_folder{:,folder_index}=velocity_state_stim;
    %%
    if plot_each_folder==1
        data_for_specific_folder_all_intensity=behavior_state_stim_folder{:,folder_index};
        
        f = figure;
        f.Position = [100 100 1520 640];
        %%% determining fraction occupancy
        for st=1:size(stim_intensities,2)
            fraction_occupancy_forward_state_folder=[];
            fraction_occupancy_turn_state_folder=[];
            fraction_occupancy_reverse_state_folder=[];
            for mn=1:size(data_for_specific_folder_all_intensity{1,1},2)
                fraction_occupancy_forward_state_folder(1,mn)=nansum(data_for_specific_folder_all_intensity{st,1}(:,mn)==1)./sum(~isnan(data_for_specific_folder_all_intensity{st,1}(:,mn)));
                fraction_occupancy_turn_state_folder(1,mn)=nansum(data_for_specific_folder_all_intensity{st,1}(:,mn)==2)./sum(~isnan(data_for_specific_folder_all_intensity{st,1}(:,mn)));
                fraction_occupancy_reverse_state_folder(1,mn)=nansum(data_for_specific_folder_all_intensity{st,1}(:,mn)==3)./sum(~isnan(data_for_specific_folder_all_intensity{st,1}(:,mn))); 
            end
            
            subplot(1,3,st)
            plot([(-1*time_window_before/fps):1/fps:time_window_after/fps], fraction_occupancy_forward_state_folder,'LineWidth',2,'Color',behavior_colors(1,:),'DisplayName',behavior_names{1})
            hold on
            plot([(-1*time_window_before/fps):1/fps:time_window_after/fps],fraction_occupancy_turn_state_folder,'LineWidth',2,'Color',behavior_colors(2,:),'DisplayName',behavior_names{2})
            hold on
            plot([(-1*time_window_before/fps):1/fps:time_window_after/fps],fraction_occupancy_reverse_state_folder,'LineWidth',2,'Color',behavior_colors(3,:),'DisplayName',behavior_names{3})
            hold on
            light_pulse = area([0 current_rails_dur/fps], [1 1],'facecolor',[1 0 0],'FaceAlpha',0.1,'LineStyle','none');
            hold on; 
            xlabel('Time (sec)')
            ylabel('Behavioral Ratio')
            ax = gca; ax.FontSize = 18;
            title(['n= ' num2str(size(data_for_specific_folder_all_intensity{st,1},1)) ', ' 'Stim intensity=' num2str(stim_intensities(1,st)) '\muW/mm^2'],'FontSize',14,'FontWeight','normal')
            ylim([0 1])
            yticks(0:0.2:1)
            set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            legend('show','Location','best');
            sgtitle(['Folder number=' num2str(folder_index) ', ' 'Total number of tracks=' num2str(size(current_tracks,2))],'fontsize',18)
        end
    end
    
end

end

%% combining the results from all the folders
behavior_state_stim_folder=behavior_state_stim_folder(~cellfun('isempty',behavior_state_stim_folder));
velocity_state_stim_folder=velocity_state_stim_folder(~cellfun('isempty',velocity_state_stim_folder));

compiled_behavior_data=cell(length(stim_intensities),1);
for frow=1:size(behavior_state_stim_folder{1,1},1) 
    dummy_data=[];
    for fcolumn=1:size(behavior_state_stim_folder,2)
        dummy_folder_data=behavior_state_stim_folder{1,fcolumn};
        dummy_data=[dummy_data;dummy_folder_data{frow,1}];
    end
    compiled_behavior_data{frow,1}=dummy_data;
end

compiled_velocity_data=cell(length(stim_intensities),1);
for frow=1:size(velocity_state_stim_folder{1,1},1) 
    dummy_data_velocity=[];
    for fcolumn=1:size(velocity_state_stim_folder,2)
        dummy_folder_data_velocity=velocity_state_stim_folder{1,fcolumn};
        dummy_data_velocity=[dummy_data_velocity;dummy_folder_data_velocity{frow,1}];
    end
    compiled_velocity_data{frow,1}=dummy_data_velocity;
end
%%%% determining fraction occupancy for all the data
%%
final_data=cell(3,3);
for mn=1:size(stim_intensities,2)
    fraction_occupancy_forward_state=[];
    fraction_occupancy_turn_state=[];
    fraction_occupancy_reverse_state=[];
    data_for_specific_intensity=compiled_behavior_data{mn,:};
    for pq=1:size(data_for_specific_intensity,2)
        fraction_occupancy_forward_state(1,pq)=nansum(data_for_specific_intensity(:,pq)==1)./sum(~isnan(data_for_specific_intensity(:,pq)));
        fraction_occupancy_turn_state(1,pq)=nansum(data_for_specific_intensity(:,pq)==2)./sum(~isnan(data_for_specific_intensity(:,pq)));
        fraction_occupancy_reverse_state(1,pq)=nansum(data_for_specific_intensity(:,pq)==3)./sum(~isnan(data_for_specific_intensity(:,pq))); 
    end
    
    final_data{mn,1}=fraction_occupancy_forward_state;
    final_data{mn,2}=fraction_occupancy_turn_state;
    final_data{mn,3}=fraction_occupancy_reverse_state;
    
end
%%
%%%% plotting the results of all the behaviors together
if plot_all_behavior_together==1
    for ijk=1:size(stim_intensities,2)
        figure;
        for beh=1:size(behavior_names,2)
            plot([(-1*time_window_before/fps):1/fps:time_window_after/fps],final_data{ijk,beh},'LineWidth',2,'LineWidth',2,'Color',behavior_colors(beh,:),'DisplayName',behavior_names{beh})
            hold on;
        end
        hold on
        light_pulse = area([0 current_rails_dur/fps], [1 1],'facecolor',[1 0 0],'FaceAlpha',0.1,'LineStyle','none');
        hold on; 
        xlabel('Time (sec)')
        ylabel('Fraction of worms')
        ax = gca; ax.FontSize = 18;
        title(['n= ' num2str(size(compiled_behavior_data{ijk,1},1)) ', ' 'Stim intensity=' num2str(stim_intensities(1,ijk)) '\muW/mm^2'],'FontSize',14,'FontWeight','normal')
        ylim([0 1])
        yticks(0:0.2:1)
        set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend('show','Location','northwest');
    end
end
%%
%%%% plotting the results of all the intensities together
if plot_all_intensities_together==1
    for jk=1:size(stim_intensities,2)
        figure;
        for beh=1:size(behavior_names,2)
            plot([(-1*time_window_before/fps):1/fps:time_window_after/fps],final_data{beh,jk},'LineWidth',2,'DisplayName',['n= ' num2str(size(compiled_behavior_data{beh,1},1)) ', ' num2str(stim_intensities(1,beh)), '\muW/mm^2 '])
            hold on;
        end

        hold on
        light_pulse = area([0 current_rails_dur/fps], [1 1],'facecolor',[1 0 0],'FaceAlpha',0.1,'LineStyle','none');
        hold on; 
        xlabel('Time (sec)')
        ylabel('Fraction of worms')
        ax = gca; ax.FontSize = 16;
        title(['Beh: ' behavior_names{jk}],'FontSize',14,'FontWeight','normal');

        ylim([0 1])
        yticks(0:0.2:1)
        set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend('show','Location','northwest');

    end
end

%%
%%%% to count for the number of reversals
count_number_of_reversal_stim=cell(3,1);
for intensity_index=1:size(stim_intensities,2)
    count_number_of_reversal_folder=[];
    for f_index=1:size(velocity_state_stim_folder,2)
        count_number_of_reversal_trial=[];
        vel_data_of_interest=velocity_state_stim_folder{1,f_index};
        dummy_vel_data_of_interest=vel_data_of_interest{intensity_index,1};
        for mn=1:size(dummy_vel_data_of_interest,1)
            dummy_vel_data_of_interest_track=dummy_vel_data_of_interest(mn,:);
            if nansum(dummy_vel_data_of_interest_track(time_window_before+1:(time_window_before+(stim_duration*fps)))<=-0.1)>7
                dummy_number_of_reversal_folder(:,1)=f_index;
                dummy_number_of_reversal_folder(:,2)=intensity_index;
                dummy_number_of_reversal_folder(:,3)=mn;
                dummy_number_of_reversal_folder(:,4)=1;
            else
                dummy_number_of_reversal_folder(:,1)=f_index;
                dummy_number_of_reversal_folder(:,2)=intensity_index;
                dummy_number_of_reversal_folder(:,3)=mn;
                dummy_number_of_reversal_folder(:,4)=0;
            end

            count_number_of_reversal_trial=[count_number_of_reversal_trial
                dummy_number_of_reversal_folder];
        end
        count_number_of_reversal_folder=[count_number_of_reversal_folder;count_number_of_reversal_trial];
    end
    count_number_of_reversal_stim{intensity_index,1}=count_number_of_reversal_folder;
end
% clc
nansum(count_number_of_reversal_stim{1, 1}(:,4))/size(count_number_of_reversal_stim{1, 1}(:,4),1);
nansum(count_number_of_reversal_stim{2, 1}(:,4))/size(count_number_of_reversal_stim{2, 1}(:,4),1);
nansum(count_number_of_reversal_stim{3, 1}(:,4))/size(count_number_of_reversal_stim{3, 1}(:,4),1);

%% Plotting barplots
if plot_barplots==1
    for i=1:size(stim_intensities,2)
        [mean_data(:,i),err_high(:,i),err_low(:,i)]=bootstrap_mean_and_ci(10000,0.05,count_number_of_reversal_stim{i, 1}(:,4));
        x_label_names{1,i}=num2str(stim_intensities(i));
    end

    x = [1:1:size(stim_intensities,2)];
    figure1=figure;
    b=bar(x,mean_data,'FaceColor','flat');                
    hold on
    er = errorbar(x,mean_data,err_low,err_high);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth=1.5;
    b.CData(1,:) = [0 0.4470 0.7410];
    b.CData(2,:) = [0.8500 0.3250 0.0980];
    b.CData(3,:) = [0.9290 0.6940 0.1250];

    set(gca,'xticklabel',x_label_names,'FontSize',14)
    % yticks([0:0.1:0.35])
    % ylim([0 0.25])
    ylabel('Fraction reversing')
    xlabel('Opto stim intensity (\muW/mm^2)')
    box off
    hold off
end
%%
%%% plotting heatmaps
% % color_matrix=[1 0 0; 0.2 1 0.2; 0 0 1];
color_matrix=behavior_colors;
% % color_matrix(1,:)=1;
if plot_behavior_heatmaps==1
    f = figure;
    f.Position = [100 100 940 420];
    for ij=1:size(compiled_behavior_data,1)
        h3=subplot(1,3,ij);
        hAxes = gca;
        clims = [0. 3.1];
        dummy_data_for_heatmap=compiled_behavior_data{ij,1};
        rng(11,'twister')
        random_samples=randsample(size(dummy_data_for_heatmap,1),60);
        imagesc(hAxes,dummy_data_for_heatmap(random_samples,:),clims)
        colormap( hAxes , color_matrix)
        xticks([1:140:281])
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-10, 10, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title([num2str(stim_intensities(ij)) '\muW/mm^2'])
        ylabel('Worm #')
        xlabel('Time (s)')
        ax = gca; ax.FontSize = 14;
        if ij ~= 1
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
           set(gca,'ylabel',[])
           set(gca,'xlabel',[])
        end
        if ij==3
            cb=colorbar;
            originalSize2 = get(gca, 'Position');
            set(h3, 'Position', originalSize2); 
        end
    end
end

%% plotting velocity plots

if plot_velocity_traces==1
    f = figure;
    f.Position = [100 100 940 420];
    time_axis=[(-1*time_window_before/fps):1/fps:time_window_after/fps];
    n_bins = 40;
    edges = linspace(-0.4,0.4,n_bins);
    velocity_density = zeros(numel(edges)-1, numel(time_axis));
    
    for ab=1:size(compiled_velocity_data,1)
       dummy_velocity_matrix=compiled_velocity_data{ab,1};
       mean_velocity=nanmean(dummy_velocity_matrix,1);
       h3=subplot(1,3,ab);
       
       for time_index = 1:length(time_axis)
           velocity_density(:,time_index) = histcounts(dummy_velocity_matrix(:,time_index), edges,'Normalization','probability');
       end
        
       hold on
       imagesc(time_axis, edges, velocity_density);
       hold on;
       plot(time_axis,mean_velocity,'-k','LineWidth',2);
       hold on;
       title([num2str(stim_intensities(ab)) '\muW/mm^2'])
       ylim([-0.4 0.4])
       yticks([-0.4:0.2:0.4])
       ylabel('Velocity (mm/sec)')
       xlabel('Time (s)')
       ax = gca; ax.FontSize = 14;
       if ab ~= 1
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            set(gca,'ylabel',[])
            set(gca,'xlabel',[])
       end
       if ab==3
            cb=colorbar;
            originalSize3 = get(gca, 'Position');
            set(h3, 'Position', originalSize3); 
       end
    end
end

%% plot velocity traces with imagesc

if plot_velocity_heatmaps==1
    
    velocity_based_behavior_colors=[0,0,1; 1,0.67,0.67; 1,0.33,0.33; 1,0,0];

    color_matrix=velocity_based_behavior_colors;
    color_matrix(2,:)=1;
    color_matrix(3,:)=1;

    f = figure;
    f.Position = [100 100 1140 520];
    for ij=1:size(compiled_velocity_data,1)
        h3=subplot(1,3,ij);
        hAxes = gca;
        clims = [-0.2 0.2];
        dummy_data_for_velocity_heatmap=compiled_velocity_data{ij,1};
        rng(1,'twister')
        random_samples=randsample(size(dummy_data_for_velocity_heatmap,1),40);
        imagesc(hAxes,dummy_data_for_velocity_heatmap(random_samples,:),clims)
        colormap( hAxes , color_matrix)
        xticks([1:140:281])
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-10, 10, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title([num2str(stim_intensities(ij)) '\muW/mm^2'])
        ylabel('Worm #')
        xlabel('Time (s)')
        ax = gca; ax.FontSize = 14;
        if ij ~= 1
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
           set(gca,'ylabel',[])
           set(gca,'xlabel',[])
        end
        if ij==3
            cb=colorbar;
            originalSize2 = get(gca, 'Position');
            set(h3, 'Position', originalSize2); 
        end

    end
end

