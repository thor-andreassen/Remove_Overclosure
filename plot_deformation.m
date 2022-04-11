%% clearing
clear
close all
clc

%% load data
load('muscle_geom_orig.mat');

%% Determine Max value
geom1_dif=geom1.vertices_history{end}-geom1.vertices_history{1};
geom1_mag=zeros(size(geom1_dif,1),1);
for count_node=1:size(geom1_dif,1)
    geom1_mag(count_node)=norm(geom1_dif(count_node,:));
end
max_deform_1=max(geom1_mag);

geom2_dif=geom2.vertices_history{end}-geom2.vertices_history{1};
geom2_mag=zeros(size(geom2_dif,1),1);
for count_node=1:size(geom2_dif,1)
    geom2_mag(count_node)=norm(geom2_dif(count_node,:));
end
max_deform_2=max(geom2_mag);

max_deform=max([max_deform_1,max_deform_2]);


%% create figure
figure('units','normalized','outerposition',[0 0 1 1])
num_frames=20;
t_vals=linspace(0,1,num_frames);

p1=patch('Faces',geom1.faces,'Vertices',geom1.vertices_history{1},'FaceColor','interp','EdgeAlpha',.3,'FaceVertexCData',geom1_mag);
colormap jet
caxis([0,max_deform]);
xlim([550,750]);
ylim([300,500]);
axis square
set(gca,'visible','off')
set(gcf,'color','w');

figure('units','normalized','outerposition',[0 0 1 1])
p2=patch('Faces',geom2.faces,'Vertices',geom2.vertices_history{1},'FaceColor','interp','EdgeAlpha',.3,'FaceVertexCData',geom2_mag);
colormap jet
caxis([0,max_deform]);
xlim([550,750]);
ylim([300,500]);
axis square
set(gca,'visible','off')
set(gcf,'color','w');


vid = VideoWriter('deform_data.gif');
open(vid);

vid2 = VideoWriter('deform_data2.gif');
open(vid2);

for count_frame=1:(numel(geom1.vertices_history)-1)
    for count_interp=t_vals
        figure(1)
        start_data_1=geom1.vertices_history{count_frame};
        end_data_1=geom1.vertices_history{count_frame+1};
        current_data_1=interp1([0,1],[reshape(start_data_1,[],1),reshape(end_data_1,[],1)]',count_interp);
        current_vertex1=reshape(current_data_1,[],3);
        p1.Vertices=current_vertex1;
        
        geom1_dif=current_vertex1-geom1.vertices_history{end};
        geom1_mag=zeros(size(geom1_dif,1),1);
        for count_node=1:size(geom1_dif,1)
            geom1_mag(count_node)=norm(geom1_dif(count_node,:));
        end
        p1.CData=geom1_mag;
        frame=getframe(gcf);
        writeVideo(vid,frame);
        
        
        figure(2)
        start_data_2=geom2.vertices_history{count_frame};
        end_data_2=geom2.vertices_history{count_frame+1};
        current_data_2=interp1([0,1],[reshape(start_data_2,[],1),reshape(end_data_2,[],1)]',count_interp);
        current_vertex2=reshape(current_data_2,[],3);
        p2.Vertices=current_vertex2;
        
        geom2_dif=current_vertex2-geom2.vertices_history{end};
        geom2_mag=zeros(size(geom2_dif,1),1);
        for count_node=1:size(geom2_dif,1)
            geom2_mag(count_node)=norm(geom2_dif(count_node,:));
        end
        p2.CData=geom2_mag;
        frame=getframe(gcf);
        writeVideo(vid2,frame);
        
    end
    
end

close(vid);
close(vid2);

%% plot still figures
clear p1 p2
close all
f1=figure('units','normalized','outerposition',[0 0 1 1])
p1=patch('Faces',geom1.faces,'Vertices',geom1.vertices_history{1},'FaceColor','interp','EdgeAlpha',.3,'FaceColor','r');
hold on
p2=patch('Faces',geom2.faces,'Vertices',geom2.vertices_history{1},'FaceColor','interp','EdgeAlpha',.3,'FaceColor','b');
xlim([550,750]);
ylim([300,500]);
colormap jet
caxis([0,max_deform]);
set(gca,'visible','off')
set(gcf,'color','w');
saveas(f1,'original_glut_data.png');

f2=figure('units','normalized','outerposition',[0 0 1 1])
p1=patch('Faces',geom1.faces,'Vertices',geom1.vertices_history{end},'FaceColor','interp','EdgeAlpha',.3,'FaceColor','r');
hold on
p2=patch('Faces',geom2.faces,'Vertices',geom2.vertices_history{end},'FaceColor','interp','EdgeAlpha',.3,'FaceColor','b');
xlim([550,750]);
ylim([300,500]);
set(gca,'visible','off')
set(gcf,'color','w');
saveas(f2,'new_glut_data.png');