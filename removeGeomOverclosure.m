%% clearing
clear
close all
clc

%% main
total_error=-Inf;
counter=1;
while total_error < -.1
%     % load stls
%     [med_geom.faces_orig, med_geom.vertices_orig]=stlRead2('glut_med_test.stl');
%     [max_geom.faces_orig, max_geom.vertices_orig]=stlRead2('glut_max_test.stl');
% 
%     save('muscle_geom_orig.mat');

    %% load mat file
    load('muscle_geom_orig.mat');

    %% Reduced Mesh
    mesh_reduction_factor=1000;
    temp=reducepatch(med_geom.faces_orig,med_geom.vertices_orig,mesh_reduction_factor);
    med_geom.faces_reduce=temp.faces;
    med_geom.vertices_reduce=temp.vertices;

    mesh_reduction_factor=1000;
    temp=reducepatch(max_geom.faces_orig,max_geom.vertices_orig,mesh_reduction_factor);
    max_geom.faces_reduce=temp.faces;
    max_geom.vertices_reduce=temp.vertices;

    %% initial plots
    med_geom_fig=figure()
    med_patch_orig=patch('Faces',med_geom.faces_orig,'Vertices',med_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    med_patch_reduce=patch('Faces',med_geom.faces_reduce,'Vertices',med_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    
    max_geom_fig=figure()
    max_patch_orig=patch('Faces',max_geom.faces_orig,'Vertices',max_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    max_patch_reduce=patch('Faces',max_geom.faces_reduce,'Vertices',max_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    
    both_geom_fig=figure()
    med_patch_orig=patch('Faces',med_geom.faces_orig,'Vertices',med_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    med_patch_reduce=patch('Faces',med_geom.faces_reduce,'Vertices',med_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    max_patch_orig=patch('Faces',max_geom.faces_orig,'Vertices',max_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    max_patch_reduce=patch('Faces',max_geom.faces_reduce,'Vertices',max_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    

    %% determine initial overclosure distances
    tri_timer=tic();
    [max_surf_to_med_point_distance, max_surf_to_med_point_points] = point2trimesh('Faces',max_geom.faces_reduce,...
        'Vertices',max_geom.vertices_reduce,'QueryPoints',med_geom.vertices_reduce,'Algorithm','parallel');

    toc(tri_timer)

    tri_timer=tic();
    [med_surf_to_max_point_distance, med_surf_to_max_point_points] = point2trimesh('Faces',med_geom.faces_reduce,...
        'Vertices',med_geom.vertices_reduce,'QueryPoints',max_geom.vertices_reduce,'Algorithm','parallel');

    toc(tri_timer)

    %% Define Gap Threshold
    desired_gap=1.5;
    relative_gap_weight=0.5; % fixed max, moving med = 0.0


    %% Determine Gaps
    max_surf_to_med_point_distance_dir=sign(max_surf_to_med_point_distance);
    max_surf_to_med_point_distance=max_surf_to_med_point_distance-desired_gap;
    max_surf_to_med_vector=max_surf_to_med_point_points;
    for count_vertex=1:length(max_surf_to_med_point_distance)
        if max_surf_to_med_point_distance(count_vertex)>0
            max_surf_to_med_point_distance(count_vertex)=0;
            multiplier=0;
        else
            multiplier=1;
        end
        max_surf_to_med_vector(count_vertex,:)=-max_surf_to_med_point_distance_dir(count_vertex)*(max_surf_to_med_point_points(count_vertex,:)...
            -med_geom.vertices_reduce(count_vertex,:))*(1-relative_gap_weight)*multiplier;
    end



    med_surf_to_max_point_distance_dir=sign(med_surf_to_max_point_distance);
    med_surf_to_max_point_distance=med_surf_to_max_point_distance-desired_gap;
    med_surf_to_max_vector=med_surf_to_max_point_points;
    for count_vertex=1:length(med_surf_to_max_point_distance)
        if med_surf_to_max_point_distance(count_vertex)>0
            med_surf_to_max_point_distance(count_vertex)=0;
            multiplier=0;
        else
            multiplier=1;
        end
        med_surf_to_max_vector(count_vertex,:)=-med_surf_to_max_point_distance_dir(count_vertex)*(med_surf_to_max_point_points(count_vertex,:)...
            -max_geom.vertices_reduce(count_vertex,:))*relative_gap_weight*multiplier;
    end

    %% Create arrow deformation plot
%     med_deform_reduce_fig=figure();
%     plot3(med_geom.vertices_reduce(:,1),med_geom.vertices_reduce(:,2),...
%         med_geom.vertices_reduce(:,3),'ro')
%     hold on
%     if norm(max_surf_to_med_vector)>0
%         arrow3(med_geom.vertices_reduce,max_surf_to_med_vector+med_geom.vertices_reduce)
%     end
%     axis equal
% 
%     max_deform_reduce_fig=figure();
%     plot3(max_geom.vertices_reduce(:,1),max_geom.vertices_reduce(:,2),...
%         max_geom.vertices_reduce(:,3),'ro')
%     hold on
%     if norm(med_surf_to_max_vector)>0
%         arrow3(max_geom.vertices_reduce,med_surf_to_max_vector+max_geom.vertices_reduce)
%     end
%     axis equal
    %% create Radial Basis Approximation
    rbf_iterations=1000;
    smoothing=10;
    rbf_timer=tic();
    med_deform_vec_rbf=newrb(med_geom.vertices_reduce',max_surf_to_med_vector',...
        1E-6,smoothing,rbf_iterations);
    toc(rbf_timer)


    rbf_timer=tic();
    max_deform_vec_rbf=newrb(max_geom.vertices_reduce',med_surf_to_max_vector',...
        1E-6,smoothing,rbf_iterations);
    toc(rbf_timer)

    %% Determine Original Deformations
    med_deform_orig_vec=sim(med_deform_vec_rbf,med_geom.vertices_orig');
    med_deform_orig_vec=med_deform_orig_vec';
    max_deform_orig_vec=sim(max_deform_vec_rbf,max_geom.vertices_orig');
    max_deform_orig_vec=max_deform_orig_vec';

    %% Plot original deformations
%     med_deform_orig_fig=figure();
%     for count_face=1:size(med_geom.faces_orig,1)
%         nodel=med_geom.faces_orig(count_face,:);
%         temp_vec=med_deform_orig_vec(nodel,:);
%         temp_vec_mag=zeros(3,1);
%         for count_vec=1:3
%             temp_vec_mag(count_vec)=norm(temp_vec(count_vec,:));
%         end
%         patch(med_geom.vertices_orig(nodel,1),med_geom.vertices_orig(nodel,2),...
%             med_geom.vertices_orig(nodel,3),temp_vec_mag);
%         hold on
%     end
%     colorbar
%     colormap jet
%     
%     
%     max_deform_orig_fig=figure();
%     for count_face=1:size(max_geom.faces_orig,1)
%         nodel=max_geom.faces_orig(count_face,:);
%         temp_vec=max_deform_orig_vec(nodel,:);
%         temp_vec_mag=zeros(3,1);
%         for count_vec=1:3
%             temp_vec_mag(count_vec)=norm(temp_vec(count_vec,:));
%         end
%         patch(max_geom.vertices_orig(nodel,1),max_geom.vertices_orig(nodel,2),...
%             max_geom.vertices_orig(nodel,3),temp_vec_mag);
%         hold on
%     end
%     colorbar
%     colormap jet


    %% Apply deformations
%     deformed_med_fig=figure()
%     plot3(med_geom.vertices_orig(:,1),med_geom.vertices_orig(:,2),...
%         med_geom.vertices_orig(:,3),'ro');
%     hold on
    med_geom.vertices_orig=med_geom.vertices_orig+med_deform_orig_vec;
%     plot3(med_geom.vertices_orig(:,1),med_geom.vertices_orig(:,2),...
%         med_geom.vertices_orig(:,3),'bo');
%     if norm(max_surf_to_med_vector)>0
%         arrow3(med_geom.vertices_reduce,max_surf_to_med_vector+med_geom.vertices_reduce)
%     end
% 
%     deformed_max_fig=figure()
%     plot3(max_geom.vertices_orig(:,1),max_geom.vertices_orig(:,2),...
%         max_geom.vertices_orig(:,3),'ro');
%     hold on
    max_geom.vertices_orig=max_geom.vertices_orig+max_deform_orig_vec;
%     plot3(max_geom.vertices_orig(:,1),max_geom.vertices_orig(:,2),...
%         max_geom.vertices_orig(:,3),'bo');
%     if norm(med_surf_to_max_vector)>0
%         arrow3(max_geom.vertices_reduce,med_surf_to_max_vector+max_geom.vertices_reduce+.0001)
%     end
    %% Display Original Min Gap
    med_error=min(max_surf_to_med_point_distance)
    max_error=min(med_surf_to_max_point_distance)
    total_error=min([med_error,max_error]);
    %% Save new stls
    stlWrite2('glut_med_test.stl',med_geom.faces_orig,med_geom.vertices_orig);
    stlWrite2('glut_max_test.stl',max_geom.faces_orig,max_geom.vertices_orig);
    
    close all
    save('muscle_geom_orig.mat');
    counter=counter+1;
end