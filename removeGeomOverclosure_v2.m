%% clearing
clear
close all
clc

%% main
total_error=-Inf;
counter=1;
while total_error < -.1
%     load stls
%     [geom2_geom.faces_orig, geom2_geom.vertices_orig]=stlRead2('glut_med_test.stl');
%     [geom1_geom.faces_orig, geom1_geom.vertices_orig]=stlRead2('glut_max_test.stl');
% 
%     save('muscle_geom_orig.mat');

    %% load mat file
    load('muscle_geom_orig.mat');

    %% Reduced Mesh
    geom2_mesh_reduction_factor=2000;
    temp=reducepatch(geom2_geom.faces_orig,geom2_geom.vertices_orig,geom2_mesh_reduction_factor);
    geom2_geom.faces_reduce=temp.faces;
    geom2_geom.vertices_reduce=temp.vertices;

    geom1_mesh_reduction_factor=2000;
    temp=reducepatch(geom1_geom.faces_orig,geom1_geom.vertices_orig,geom1_mesh_reduction_factor);
    geom1_geom.faces_reduce=temp.faces;
    geom1_geom.vertices_reduce=temp.vertices;

    %% initial plots
    geom2_geom_fig=figure()
    geom2_patch_orig=patch('Faces',geom2_geom.faces_orig,'Vertices',geom2_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    geom2_patch_reduce=patch('Faces',geom2_geom.faces_reduce,'Vertices',geom2_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    
    geom1_geom_fig=figure()
    geom1_patch_orig=patch('Faces',geom1_geom.faces_orig,'Vertices',geom1_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    geom1_patch_reduce=patch('Faces',geom1_geom.faces_reduce,'Vertices',geom1_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    
    both_geom_fig=figure()
    geom2_patch_orig=patch('Faces',geom2_geom.faces_orig,'Vertices',geom2_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    geom2_patch_reduce=patch('Faces',geom2_geom.faces_reduce,'Vertices',geom2_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    geom1_patch_orig=patch('Faces',geom1_geom.faces_orig,'Vertices',geom1_geom.vertices_orig,'FaceColor','r','FaceAlpha',.5)
    hold on
    geom1_patch_reduce=patch('Faces',geom1_geom.faces_reduce,'Vertices',geom1_geom.vertices_reduce,'FaceColor','g','FaceAlpha',.5)
    

    %% determine initial overclosure distances
    tri_timer=tic();
    [geom1_surf_to_geom2_point_distance, geom1_surf_to_geom2_point_points] = point2trimesh('Faces',geom1_geom.faces_reduce,...
        'Vertices',geom1_geom.vertices_reduce,'QueryPoints',geom2_geom.vertices_reduce,'Algorithm','parallel');

    toc(tri_timer)

    tri_timer=tic();
    [geom2_surf_to_geom1_point_distance, geom2_surf_to_geom1_point_points] = point2trimesh('Faces',geom2_geom.faces_reduce,...
        'Vertices',geom2_geom.vertices_reduce,'QueryPoints',geom1_geom.vertices_reduce,'Algorithm','parallel');

    toc(tri_timer)

    %% Define Gap Threshold
%     desired_gap=1.5;
%     relative_gap_weight=0.5% fixed geom1, moving geom2 = 0.0
    desired_gap=15;
    relative_gap_weight=0.5; % fixed geom1, moving geom2 = 0.0


    %% Determine Gaps
    geom1_surf_to_geom2_point_distance_dir=sign(geom1_surf_to_geom2_point_distance);
    geom1_surf_to_geom2_point_distance=geom1_surf_to_geom2_point_distance-desired_gap;
    geom1_surf_to_geom2_vector=geom1_surf_to_geom2_point_points;
    for count_vertex=1:length(geom1_surf_to_geom2_point_distance)
        if geom1_surf_to_geom2_point_distance(count_vertex)>0
            geom1_surf_to_geom2_point_distance(count_vertex)=0;
            multiplier=0;
        else
            multiplier=1;
        end
        geom1_surf_to_geom2_vector(count_vertex,:)=-geom1_surf_to_geom2_point_distance_dir(count_vertex)*(geom1_surf_to_geom2_point_points(count_vertex,:)...
            -geom2_geom.vertices_reduce(count_vertex,:))*(1-relative_gap_weight)*multiplier;
    end



    geom2_surf_to_geom1_point_distance_dir=sign(geom2_surf_to_geom1_point_distance);
    geom2_surf_to_geom1_point_distance=geom2_surf_to_geom1_point_distance-desired_gap;
    geom2_surf_to_geom1_vector=geom2_surf_to_geom1_point_points;
    for count_vertex=1:length(geom2_surf_to_geom1_point_distance)
        if geom2_surf_to_geom1_point_distance(count_vertex)>0
            geom2_surf_to_geom1_point_distance(count_vertex)=0;
            multiplier=0;
        else
            multiplier=1;
        end
        geom2_surf_to_geom1_vector(count_vertex,:)=-geom2_surf_to_geom1_point_distance_dir(count_vertex)*(geom2_surf_to_geom1_point_points(count_vertex,:)...
            -geom1_geom.vertices_reduce(count_vertex,:))*relative_gap_weight*multiplier;
    end

    %% Create arrow deformation plot
%     geom2_deform_reduce_fig=figure();
%     plot3(geom2_geom.vertices_reduce(:,1),geom2_geom.vertices_reduce(:,2),...
%         geom2_geom.vertices_reduce(:,3),'ro')
%     hold on
%     if norm(geom1_surf_to_geom2_vector)>0
%         arrow3(geom2_geom.vertices_reduce,geom1_surf_to_geom2_vector+geom2_geom.vertices_reduce)
%     end
%     axis equal
% 
%     geom1_deform_reduce_fig=figure();
%     plot3(geom1_geom.vertices_reduce(:,1),geom1_geom.vertices_reduce(:,2),...
%         geom1_geom.vertices_reduce(:,3),'ro')
%     hold on
%     if norm(geom2_surf_to_geom1_vector)>0
%         arrow3(geom1_geom.vertices_reduce,geom2_surf_to_geom1_vector+geom1_geom.vertices_reduce)
%     end
%     axis equal
    %% create Radial Basis Approximation
%     rbf_iterations=1000;
%     smoothing=10;
    rbf_iterations=1000;
    smoothing=1000;
    rbf_timer=tic();
%     geom2_deform_vec_rbf=newrb(geom2_geom.vertices_reduce',geom1_surf_to_geom2_vector',...
%         1E-6,smoothing,rbf_iterations);
%     toc(rbf_timer)
% 
% 
%     rbf_timer=tic();
%     geom1_deform_vec_rbf=newrb(geom1_geom.vertices_reduce',geom2_surf_to_geom1_vector',...
%         1E-6,smoothing,rbf_iterations);
%     toc(rbf_timer)


    geom2_deform_vec_rbf=newrb([geom2_geom.vertices_reduce',geom1_geom.vertices_reduce'],...
        [geom1_surf_to_geom2_vector',-geom2_surf_to_geom1_vector'],...
        1E-6,smoothing,rbf_iterations);
    toc(rbf_timer)


    geom1_deform_vec_rbf=newrb([geom1_geom.vertices_reduce',geom2_geom.vertices_reduce'],...
        [geom2_surf_to_geom1_vector',-geom1_surf_to_geom2_vector'],...
        1E-6,smoothing,rbf_iterations);
    toc(rbf_timer)


    %% Determine Original Deformations
    geom2_deform_orig_vec=sim(geom2_deform_vec_rbf,geom2_geom.vertices_orig');
    geom2_deform_orig_vec=geom2_deform_orig_vec';
    geom1_deform_orig_vec=sim(geom1_deform_vec_rbf,geom1_geom.vertices_orig');
    geom1_deform_orig_vec=geom1_deform_orig_vec';

    %% Plot original deformations
%     geom2_deform_orig_fig=figure();
%     for count_face=1:size(geom2_geom.faces_orig,1)
%         nodel=geom2_geom.faces_orig(count_face,:);
%         temp_vec=geom2_deform_orig_vec(nodel,:);
%         temp_vec_mag=zeros(3,1);
%         for count_vec=1:3
%             temp_vec_mag(count_vec)=norm(temp_vec(count_vec,:));
%         end
%         patch(geom2_geom.vertices_orig(nodel,1),geom2_geom.vertices_orig(nodel,2),...
%             geom2_geom.vertices_orig(nodel,3),temp_vec_mag);
%         hold on
%     end
%     colorbar
%     colormap jet
%     
%     
%     geom1_deform_orig_fig=figure();
%     for count_face=1:size(geom1_geom.faces_orig,1)
%         nodel=geom1_geom.faces_orig(count_face,:);
%         temp_vec=geom1_deform_orig_vec(nodel,:);
%         temp_vec_mag=zeros(3,1);
%         for count_vec=1:3
%             temp_vec_mag(count_vec)=norm(temp_vec(count_vec,:));
%         end
%         patch(geom1_geom.vertices_orig(nodel,1),geom1_geom.vertices_orig(nodel,2),...
%             geom1_geom.vertices_orig(nodel,3),temp_vec_mag);
%         hold on
%     end
%     colorbar
%     colormap jet


    %% Apply deformations
%     deforgeom2_geom2_fig=figure()
%     plot3(geom2_geom.vertices_orig(:,1),geom2_geom.vertices_orig(:,2),...
%         geom2_geom.vertices_orig(:,3),'ro');
%     hold on
    geom2_geom.vertices_orig=geom2_geom.vertices_orig+geom2_deform_orig_vec;
%     plot3(geom2_geom.vertices_orig(:,1),geom2_geom.vertices_orig(:,2),...
%         geom2_geom.vertices_orig(:,3),'bo');
%     if norm(geom1_surf_to_geom2_vector)>0
%         arrow3(geom2_geom.vertices_reduce,geom1_surf_to_geom2_vector+geom2_geom.vertices_reduce)
%     end
% 
%     deforgeom2_geom1_fig=figure()
%     plot3(geom1_geom.vertices_orig(:,1),geom1_geom.vertices_orig(:,2),...
%         geom1_geom.vertices_orig(:,3),'ro');
%     hold on
    geom1_geom.vertices_orig=geom1_geom.vertices_orig+geom1_deform_orig_vec;
%     plot3(geom1_geom.vertices_orig(:,1),geom1_geom.vertices_orig(:,2),...
%         geom1_geom.vertices_orig(:,3),'bo');
%     if norm(geom2_surf_to_geom1_vector)>0
%         arrow3(geom1_geom.vertices_reduce,geom2_surf_to_geom1_vector+geom1_geom.vertices_reduce+.0001)
%     end
    %% Display Original Min Gap
    geom2_error=min(geom1_surf_to_geom2_point_distance);
    geom1_error=min(geom2_surf_to_geom1_point_distance);
    table(geom2_error,geom1_error)
    total_error=min([geom2_error,geom1_error]);
    %% Save new stls
    stlWrite2('glut_geom2_test.stl',geom2_geom.faces_orig,geom2_geom.vertices_orig);
    stlWrite2('glut_geom1_test.stl',geom1_geom.faces_orig,geom1_geom.vertices_orig);
    
    close all
    save('muscle_geom_orig.mat');
    counter=counter+1;
end