%% clearing
clear
close all
clc

%% get list of stls

%% original test path
stl_folder='C:\Users\DU_P620\Desktop\Thors_Personal_Folder\Example_Overclosure\';
result_folder=[stl_folder,'fixed_stls\'];
stl_files=dir([stl_folder,'*.stl']);
stl_file_names={stl_files.name};

%% create list of possible intersections
%% load geoms
% 
% [nodelist1,elemlist1,elemlist1_renum]=READ_MESH_NUMS_AJC('S193761_Left_Cartilage_Femur_Hex_Mesh.inp');
% geom1.faces=elemlist1_renum(:,2:end);
% geom1.vertices=nodelist1(:,2:end);
% geom1.vertices_orig=nodelist1(:,2:end);
% geom1.elemlist=elemlist1;
% geom1.nodelist=nodelist1;


% [nodelist2,elemlist2,elemlist_renum2]=READ_MESH_NUMS_AJC('S193761_Left_Cartilage_Tibia_Medial_Hex_Mesh.inp');
% geom2.faces=elemlist_renum2(:,2:end);
% geom2.vertices=nodelist2(:,2:end);
% geom2.vertices_orig=nodelist2(:,2:end);
% geom2.elemlist=elemlist2;
% geom2.nodelist=nodelist2;
% % % 
% % % intersection_matrix=zeros(length(stl_file_names));
% % % distance_matrix=intersection_matrix;
% % % 
% save('cart_geom_orig.mat');
load('cart_geom_orig.mat');
%% overclosure parameters
% % working parameters
%         params.desired_gap=.05;
%         params.stop_tolerance=1E-5;
%         params.relative_gap_weight=0.5;
%         params.element_3d_type=[0,0];
%         params.use_parallel_loops=1;
%         params.smoothing_improve=.1;
%         params.plot_surf=0;
%         params.smoothing=0.999;
%         params.smoothing_reduction=0.9;
%         params.rbf_iterations=400;
%         params.geom1_mesh_reduction_factor=1;
%         params.geom2_mesh_reduction_factor=1;
%         params.scale_percent_factor=1.3;

        params.desired_gap=.05;
        params.stop_tolerance=1E-5;
        params.relative_gap_weight=0.5;
        params.element_3d_type=[1,1];
        params.use_parallel_loops=1;
        params.smoothing_improve=1000;
        params.plot_surf=1;
        params.smoothing=1000;
%         params.smoothing=1000;
        params.smoothing_reduction=0.999;
%         params.smoothing_reduction=0.99;
        params.rbf_iterations=4000;
        params.geom1_mesh_reduction_factor=1;
        params.geom2_mesh_reduction_factor=1;
        params.scale_percent_factor=1.3;
        


%% Main over-closure adjustment loop

    try
        [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure]=...
            removeOverclosureRBF(geom1,geom2,params);
    catch
        geom1_new=geom1;
        geom2_new=geom2;
        counter=1000;
        original_max_overclosure_1=1000;
        original_max_overclosure_2=1000;
        original_max_overclosure=1000;
    end
geom1=geom1_new;
geom2=geom2_new;
save('cart_geom_orig.mat');