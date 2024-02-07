%% clearing
clear
close all
clc

%% Base Parameters
%% load geoms
[nodelist1,elemlist1,elemlist1_renum]=READ_MESH_NUMS_AJC('S193761_Left_Cartilage_Femur_Hex_Mesh.inp');
geom1.faces=elemlist1_renum(:,2:end);
geom1.vertices=nodelist1(:,2:end);
geom1.vertices_orig=nodelist1(:,2:end);
geom1.elemlist=elemlist1;
geom1.nodelist=nodelist1;


[nodelist2,elemlist2,elemlist_renum2]=READ_MESH_NUMS_AJC('S193761_Left_Cartilage_Tibia_Lateral_Hex_Mesh.inp');
geom2.faces=elemlist_renum2(:,2:end);
geom2.vertices=nodelist2(:,2:end);
geom2.vertices_orig=nodelist2(:,2:end);
geom2.elemlist=elemlist2;
geom2.nodelist=nodelist2;

% save('cart_geom_orig.mat');
load('cart_geom_orig.mat');
%% overclosure parameters
        params.desired_gap=.1;
        params.stop_tolerance=1E-5;
        params.relative_gap_weight=0.5;
        params.element_3d_type=[1,1];
        params.use_parallel_loops=1;
        params.smooth_2D_surface=500;
        params.plot_surf=1;
        params.smoothing=10;
        params.smoothing_reduction=0.99995;
        params.geom1_mesh_reduction_factor=.001;
        params.geom2_mesh_reduction_factor=.001;
        params.scale_reduction_factor=1.005;
        params.weight_factor=100;
        params.accelerated_weight=1;
        params.check_original=1;

%% Main over-closure adjustment loop

    try
        % the following line uses the GRNN algorithm (THIS IS THE MAIN
        % FUNCTION)
        [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure,history_params]=...
            removeOverclosureGRNN(geom1,geom2,params);

        % the following line uses the conventional nodal adjustment
        % algorithm
%         [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure,history_params]=...
%             removeOverclosureNODAL(geom1,geom2,params)
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
save('cart_geom_fixed_GRNN.mat');




%% output data

params.elem_type='C3D8R';
params.nset_name='FEMUR_CART_NODES';
params.elset_name='FEMUR_CART_ELEMS';

nodes1=geom1.nodelist;
nodes1(:,2:end)=geom1.vertices;
elements1=geom1.elemlist;

writeAbaqusInput('S193761_Left_Cartilage_Femur_Hex_Mesh_no_over_validation_GRNN.inp',nodes1,elements1,params);

params.elem_type='C3D8R';
params.nset_name='TIBIA_CART_LAT_NODES';
params.elset_name='TIBIA_CART_LAT_ELEMS';

nodes1=geom2.nodelist;
nodes1(:,2:end)=geom2.vertices;
elements1=geom2.elemlist;

writeAbaqusInput('S193761_Left_Cartilage_Tibia_Lateral_Hex_Mesh_no_over_validation_GRNN.inp',nodes1,elements1,params);