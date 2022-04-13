%% clearing
clear
close all
clc

%% get list of stls

%% original test path
stl_folder='C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\visible human\overclosure corrected stls\STLs\Left\';


stl_files=dir([stl_folder,'*.stl']);
stl_file_names={stl_files.name};

%% create list of possible intersections
threshold_dist=1;
intersection_matrix=zeros(length(stl_file_names));
for counti=1:length(stl_file_names)
    counti
        [geom1_temp.faces,geom1_temp.vertices]=stlRead2([stl_folder,stl_file_names{counti}]);
        geom1=reducepatch(geom1_temp.faces,geom1_temp.vertices,1000);
        for countj=counti:length(stl_file_names)
                [geom2_temp.faces,geom2_temp.vertices]=stlRead2([stl_folder,stl_file_names{countj}]);
                geom2=reducepatch(geom2_temp.faces,geom2_temp.vertices,1000);
               [ distances1 ] = point2trimesh('Faces',geom1.faces,...
                       'Vertices',geom1.vertices,'QueryPoints',geom2.vertices,'Algorithm','parallel');
               [ distances2 ] = point2trimesh('Faces',geom2.faces,...
                       'Vertices',geom2.vertices,'QueryPoints',geom1.vertices,'Algorithm','parallel');
               min1=min(distances1);
               min2=min(distances2);
               min_total=min([min1,min2]);
               if min_total<= threshold_dist
                      intersection_matrix(counti,countj)=1;
               end
        end
end

%% set diagonal to zero
for counti=1:length(stl_file_names)
        intersection_matrix(counti,counti)=0;
end

%% create list of possible overclsoures
    geom_name_list=[];
    counter=1;
    for counti=1:length(stl_file_names)
        for countj=counti:length(stl_file_names)
            if intersection_matrix(counti,countj)==1
                geom_name_list(counter).geom1=stl_file_names{counti};
                geom_name_list(counter).geom2=stl_file_names{countj};
                counter=counter+1;
            end
        end
    end
    
    
    %% create overclosure job setup
    overclosure_job_list=geom_name_list;
    temp_val=size(geom_name_list);
    num_possible_over=temp_val(2);
    for count_pair=1:num_possible_over
        if testCharPresentInChar(geom_name_list(count_pair).geom1,'Bone',0)
            if testCharPresentInChar(geom_name_list(count_pair).geom2,'Bone',0)
                overclosure_job_list(count_pair).weight=0.5;
            else
                overclosure_job_list(count_pair).weight=1;
            end
        elseif testCharPresentInChar(geom_name_list(count_pair).geom2,'Bone',0)
            overclosure_job_list(count_pair).weight=0;
        else
            overclosure_job_list(count_pair).weight=0.5;
        end
        
        if testCharPresentInChar(geom_name_list(count_pair).geom1,'Bone',0) ||...
                testCharPresentInChar(geom_name_list(count_pair).geom2,'Bone',0)
            overclosure_job_list(count_pair).priority=4;
        elseif testCharPresentInChar(geom_name_list(count_pair).geom1,'Cartilage',0) ||...
                testCharPresentInChar(geom_name_list(count_pair).geom2,'Cartilage',0)
            overclosure_job_list(count_pair).priority=3;
        elseif testCharPresentInChar(geom_name_list(count_pair).geom1,'Ligament',0) ||...
                testCharPresentInChar(geom_name_list(count_pair).geom2,'Ligament',0)
            overclosure_job_list(count_pair).priority=2;
        else
            overclosure_job_list(count_pair).priority=1;
        end
    end
    
T = struct2table(overclosure_job_list);
sortedT = sortrows(T, 'priority');
overclosure_job_list = table2struct(sortedT);
    
%% overclosure parameters
        params.desired_gap=.05;
        params.relative_gap_weight=0.5;
        params.element_3d_type=[0,0];
        params.use_parallel_loops=1;
        params.smoothing_improve=.1;
        params.plot_surf=0;
        params.smoothing=0.999;
        params.rbf_iterations=400;
        params.geom1_mesh_reduction_factor=1;
        params.geom2_mesh_reduction_factor=1;
        params.scale_percent_factor=1.3;

%% Main over-closure adjustment loop
result_folder=[stl_folder,'fixed_stls\'];
w=waitbar(0,'Adjusting Over-closures');
for count_pair=72:num_possible_over
    tic
    geom1_name=overclosure_job_list(count_pair).geom1
    geom2_name=overclosure_job_list(count_pair).geom2
    if isfile([result_folder,geom1_name])
        current_path_geom1=[result_folder,geom1_name];
    else
        current_path_geom1=[stl_folder,geom1_name];
    end
    
    if isfile([result_folder,geom2_name])
        current_path_geom2=[result_folder,geom2_name];
    else
        current_path_geom2=[stl_folder,geom2_name];
    end
    
    [geom1.faces,geom1.vertices]=stlRead2(current_path_geom1);
    [geom2.faces,geom2.vertices]=stlRead2(current_path_geom2);
    
    params.relative_gap_weight=overclosure_job_list(count_pair).weight;
    
    try
        [geom1_new,geom2_new,counter,original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure]=...
            removeOverclosureRBF(geom1,geom2,params);
    catch
        geom1_new=geom1;
        geom2_new=geom2;
        counter=1000;
        original_max_overclosure_1=1000;
        original_max_overclosure_2=1000;
        original_max_overclosure_3=1000;
    end
    toc
    close all
    current_cells={geom1_name,geom2_name,counter,...
        original_max_overclosure_1,original_max_overclosure_2,original_max_overclosure};
    
    stlwrite([result_folder,geom1_name],geom1_new);
    stlwrite([result_folder,geom2_name],geom2_new);
    
    
    if count_pair==1
        writecell(current_cells,[result_folder,'log_file.xlsx']);
    else
        writecell(current_cells,[result_folder,'log_file.xlsx'],'WriteMode','append');
    end
    waitbar(count_pair/num_possible_over,w,'Adjusting Over-closures');
end