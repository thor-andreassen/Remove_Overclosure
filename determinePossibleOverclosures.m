%% clearing
clear
close all
clc

%% get list of stls
stl_files=dir('*.stl');
stl_file_names={stl_files.name};

%% create list of possible intersections
threshold_dist=10;


intersection_test=zeros(length(stl_file_names));
for counti=1:length(stl_file_names)
        geom1_temp=stlread(stl_file_names{counti});
        geom1=reducepatch(geom1_temp.ConnectivityList,geom1_temp.Points,1000);
        for countj=counti:length(stl_file_names)
                geom2_temp=stlread(stl_file_names{countj});
                geom2=reducepatch(geom2_temp.ConnectivityList,geom2_temp.Points,1000);
               [ distances1 ] = point2trimesh('Faces',geom1.faces,...
                       'Vertices',geom1.vertices,'QueryPoints',geom2.vertices,'Algorithm','parallel');
               [ distances2 ] = point2trimesh('Faces',geom2.faces,...
                       'Vertices',geom2.vertices,'QueryPoints',geom1.vertices,'Algorithm','parallel');
               min1=min(distances1);
               min2=min(distances2);
               min_total=min([min1,min2]);
               if min_total<= threshold_dist
                      intersection_test(counti,countj)=1;
               end
        end
end

%% set diagonal to zero
for counti=1:length(stl_file_names)
        intersection_test(counti,counti)=0;
end