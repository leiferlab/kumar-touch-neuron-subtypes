function [] = plot_worm_frame_SK(I, center_lines, centerline_color,plotting_index,in_track_index,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode,fps)
    
    debugimage=plotting_index;
    I=imadjust(I);
     
    if nargin > 8
        %debugimage inputted
        imshow(I + debugimage, [], 'InitialMagnification', 300, 'Border','tight');
    else
        imshow(I, [], 'InitialMagnification', 300, 'Border','tight');
    end
    
    %%% path
    smoothed_path=[smoothed_path_x smoothed_path_y];
    
    if isempty(centerline_color)
        centerline_color = [1, 1, 0];
    end
    
    centerline_correction=fliplr(center_lines(10,:));
    %%%% plot the previous few seconds of path or else plot from the beginning
    if in_track_index>path_history
        path_t = smoothed_path(in_track_index-path_history:in_track_index,:) - smoothed_path(in_track_index,:)+centerline_correction;
    else
        path_t = smoothed_path(1:in_track_index,:) - smoothed_path(in_track_index,:)+centerline_correction;
    end

    if nargin > 1
        hold on
        plot(center_lines(:,2), center_lines(:,1), '-', 'Color', centerline_color, 'LineWidth',4)
        plot(center_lines(1,2), center_lines(1,1), '.', 'Color', centerline_color, 'markersize',40)
        hold off;
    end
    hold on;
    plot(path_t(:,1),path_t(:,2),'.-g','LineWidth',4);
    hold on;
    
    %%% drawing the points separated by a set number
    plot(path_t(1:fps:end,1),path_t(1:fps:end,2),'.', 'MarkerSize',18,'color',[1 1 0],'LineWidth',2);
    hold on;
    if demo_mode~=1
    %score
    text(10, 10, num2str(in_track_index), 'Color', 'y','FontSize',14);
% %     text(60, 10, num2str(round(abs(stimulus_delivered))), 'Color', 'y','FontSize',14)
% %     %eccentricity
% %     text(10, 60, num2str(current_ellipse_ratio,'%4.2f'), 'Color', 'y','FontSize',14);
% %     %current velocity
% %     text(60, 60, num2str(current_velocity,'%4.2f'), 'Color', 'y','FontSize',14);
    %%% stim intensity
    text(120,10, [num2str(round(abs(stimulus_delivered))), ' uW/mm^{2}'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'color', 'y', 'fontsize', 14);
    %%% plotting stim circle only when there is a stim
    if stimulus_delivered>=10
        plot(100, 25, 'o', 'MarkerSize', 30, 'MarkerEdgeColor','none', 'MarkerFaceColor',[min(max(stimulus_delivered/200,0),1) 0 0])
    end
    
    hold off;
    end
end