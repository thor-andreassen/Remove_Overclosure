%% clearing
clear
close all
clc

%% load geoms
geom1_filename='glut_max_test.stl';
[geom1.faces,geom1.vertices]=stlRead2(geom1_filename);
geom1.vertices_orig=geom1.vertices;

% [nodelist,elemlist,elemlist_renum]=READ_MESH_NUMS_AJC('CART1-FEMUR-S192803L_2.inp');
% [nodelist,elemlist,elemlist_renum]=READ_MESH_NUMS_AJC('hex_cylinder.inp');
% geom1.faces=elemlist_renum(:,2:end);
% geom1.vertices=nodelist(:,2:end);
% geom1.vertices_orig=nodelist(:,2:end);

geom2_filename='glut_med_test.stl';
[geom2.faces,geom2.vertices]=stlRead2(geom2_filename);
geom2.vertices_orig=geom2.vertices;
% [nodelist,elemlist,elemlist_renum]=READ_MESH_NUMS_AJC('BONE1-FEMUR-S192803L.inp');
% [nodelist,elemlist,elemlist_renum]=READ_MESH_NUMS_AJC('tri_cone.inp');
% geom2.faces=elemlist_renum(:,2:end);
% geom2.vertices=nodelist(:,2:end);
% geom2.vertices_orig=nodelist(:,2:end);
% save('muscle_geom_orig.mat');

% load('muscle_geom_orig.mat');

%% create overclosure
desired_over=-5;
min_over=Inf;
count=1;
max_iter=100;
while count<max_iter && min_over-desired_over>.0001
    count
   [geom1_surf_to_geom2_point_distance, geom1_surf_to_geom2_point_points] = point2trimesh('Faces',geom1.faces,...
            'Vertices',geom1.vertices,'QueryPoints',geom2.vertices,'Algorithm','parallel');
   [min_over_1,index]=min(geom1_surf_to_geom2_point_distance);
   min_over_dir_1=geom1_surf_to_geom2_point_points(index,:)-geom2.vertices(index,:);
   
   [geom2_surf_to_geom1_point_distance, geom2_surf_to_geom1_point_points] = point2trimesh('Faces',geom2.faces,...
            'Vertices',geom2.vertices,'QueryPoints',geom1.vertices,'Algorithm','parallel');
   [min_over_2,index]=min(geom2_surf_to_geom1_point_distance);
   min_over_dir_2=geom2_surf_to_geom1_point_points(index,:)-geom1.vertices(index,:);
   if min_over_1<min_over_2
       min_over=min_over_1;
       u_vec=min_over_dir_1/norm(min_over_dir_1);
       geom2.vertices=geom2.vertices+(abs(min_over-desired_over)*u_vec);
   else
       min_over=min_over_2;
       u_vec=min_over_dir_2/norm(min_over_dir_2);
       geom1.vertices=geom1.vertices+(abs(min_over-desired_over)*u_vec);
   end
    count=count+1;
end

%% export geom
[~,geom1_name]=fileparts(geom1_filename);
[~,geom2_name]=fileparts(geom2_filename);


stlwrite([char(geom1_name),'_5d0.stl'],geom1);
stlwrite([char(geom2_name),'_5d0.stl'],geom2);