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
    
%% remove geometry over-closure
% [geom1_new,geom2_new]=removeOverclosureRBF(geom1,geom2);
