%%% To plot individual worm videos for tap assay setup
%%% April 5, 2022

clear
clc
close all

main_folder=('/projects/LEIFER/Sandeep/Data/20220328_array_light_mix_tap_duration_details_Sandeep_CGZ195/Data20220328_183810'); %% main folder path
load(fullfile(main_folder,'workspace_data_20220328_183810_3sec.mat')); %% save the workspace data inside the main folder

main_folder=('/projects/LEIFER/Sandeep/Data/20220328_array_light_mix_tap_duration_details_Sandeep_CGZ195/Data20220328_183810'); %% main folder path

%%%%%%%%%%%%%%%%%%%%%%%%%%%% user inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_video=1;
output_video_frame_rate=20;
fps=14;
frames_per_plot_time = round(fps/output_video_frame_rate);
plotting_index = 1;
path_history=10*fps;
demo_mode=11;

addpath(fullfile(main_folder,'/individual_worm_imgs'));
worm_frame_matrix=locs_track_index_array;

figure1=figure('Position', [40 40 820 500]);

for num_iteration=1:size(worm_frame_matrix,1) 
% for num_iteration=2 %%%%% for debugging

worm_of_interest=worm_frame_matrix(num_iteration,1);
sprintf('Analyzing worm %d out of %d worms', num_iteration, size(worm_frame_matrix,1))

str1 = "worm_";
str2 = num2str(worm_of_interest);
str3 = ".mat";

worm_image_data = append(str1,str2,str3);
load(worm_image_data)

%%%%%%%%%%%% plotting individual worm videos %%%%%%%%%%%%%%
initial_frame=worm_frame_matrix(num_iteration,2);

if save_video==1
    if ~exist([main_folder, filesep, 'individual_worm_videos/rails_200uW_3sec_vel_subplot'], 'dir')
        mkdir([main_folder, filesep, 'individual_worm_videos/rails_200uW_3sec_vel_subplot'])
    end
    video_file_name = fullfile([main_folder, filesep, 'individual_worm_videos/rails_200uW_3sec_vel_subplot', filesep, 'worm_', num2str(worm_of_interest), '_frame_',num2str(initial_frame),'.avi']);   
    outputVideo = VideoWriter(fullfile([main_folder, filesep, 'individual_worm_videos/rails_200uW_3sec_vel_subplot', filesep, 'worm_', num2str(worm_of_interest), '_frame_',num2str(initial_frame)]),'Motion JPEG AVI');
    outputVideo.FrameRate = output_video_frame_rate;
    
    open(outputVideo);
end
        
for in_track_index=(initial_frame-140):1:(initial_frame+140)
    
    if (initial_frame-140)<=0 %%%% skip the iteration if the stim happens at the beginning of the video
        continue
    end
    
    if (initial_frame+140)>size(current_tracks(worm_of_interest).Frames,2) %%%% skip the iteration if the stim happens at the end of the video
        continue
    end

    smoothed_path_x=smooth(current_tracks(worm_of_interest).Path(:,1),10);  %%% smoothing path with window 10
    smoothed_path_y=smooth(current_tracks(worm_of_interest).Path(:,2),10);
    stimulus_delivered_all=current_tracks(worm_of_interest).LEDVoltages.*folder_paras.VoltageToPower;
    stimulus_delivered=stimulus_delivered_all(in_track_index);
    velocity_on_current_frame=current_tracks(worm_of_interest).Velocity(in_track_index);
    all_velocity=current_tracks(worm_of_interest).Velocity((initial_frame-140):(initial_frame+140),:);
    stim_color=[1 0 0];
    
    I = squeeze(worm_images(:,:,in_track_index));
    current_behavior = current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index);
    if current_behavior < 1 || current_behavior > length(behavior_colors)
        current_behavior_color = [1 1 1];
    else
        current_behavior_color = behavior_colors(current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index),:);
    end
                 
    %%%%% plot videos frames with subplots  
    plot_worm_frame_with_vel_subplot_SK(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
                plotting_index,in_track_index,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode,fps,velocity_on_current_frame,all_velocity,initial_frame,current_rails_dur,stim_color); 
    
% %     %%%%% plot videos frames 
% %     plot_worm_frame_SK(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% %                 plotting_index,in_track_index,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode,fps,velocity_on_current_frame); 
    
    hold on;
    hold off;
    if save_video==1
    writeVideo(outputVideo, getframe(gcf));
    end
end

if save_video==1
    close(outputVideo);
    video_transcode(video_file_name);
end
end

