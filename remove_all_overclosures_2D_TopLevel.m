% Created by Thor E. Andreassen, PhD
% Last Edited 2/8/2024
% This code shows an example of removing overclosures between 2 meshes
% with 3D hexahedral elements. The geometries are taken from
% a set of meshes included as part of other work.
%% clearing
clear
close all
clc

%% Base Parameters
%% load geoms
% [geom1.elems,geom1.vertices]=stlRead2('VHFL_Muscle_TibialisPosterior.stl');
% geom1.elemlist=geom1.elems;
% 
% 
% [geom2.elems,geom2.vertices]=stlRead2('VHFL_Muscle_TibialisAnterior.stl');
% geom2.elemlist=geom2.elems;

% save('muscle_orig.mat');
load('muscle_orig.mat');
%% overclosure parameters
params.desired_gap=.1;
params.stop_tolerance=1E-5;
params.relative_gap_weight=0.5;
params.element_3d_type=[0,0];
params.use_parallel_loops=1;
params.smooth_2D_surface=500;
params.plot_surf=1;
params.smoothing=100;
params.smoothing_reduction=0.99995;
params.geom1_mesh_reduction_factor=.01;
params.geom2_mesh_reduction_factor=.01;
params.scale_reduction_factor=1.005;
params.weight_factor=100;
params.accelerated_weight=1;
params.check_original=0;

%% Main over-closure adjustment loop

% try
    % the following line uses the GRNN algorithm (THIS IS THE MAIN
    % FUNCTION)
%     [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure,history_params]=...
%         removeOverclosureGRNN(geom1,geom2,params);

    % the following line uses the conventional nodal adjustment
    % algorithm
            [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure,history_params]=...
                removeOverclosureNODAL(geom1,geom2,params)
% catch
%     geom1_new=geom1;
%     geom2_new=geom2;
%     counter=1000;
%     original_max_overclosure_1=1000;
%     original_max_overclosure_2=1000;
%     original_max_overclosure=1000;
% end
geom1=geom1_new;
geom2=geom2_new;
save('muscle_fixed_NODAL.mat');
