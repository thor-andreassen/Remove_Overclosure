%% clear
clear
close all
clc


%% load data
load('muscle_geom_orig.mat');

geom_temp.faces=[geom1.faces;geom2.faces+size(geom1.vertices,1)];
geom_temp.vertices=[geom1.vertices;geom2.vertices];
geom_temp.vertices_orig=[geom1.vertices_orig;geom2.vertices_orig];

%% determine mesh qualtiy
aspects=zeros(size(geom_temp.faces,1),1);
skewness=zeros(size(geom_temp.faces,1),1);
for count_face=1:size(geom_temp.faces,1)
    nodel=geom_temp.faces(count_face,:);
    face_nodes=geom_temp.vertices(nodel,:);
    [skewness(count_face),aspects(count_face)]=getMeshQuality2(face_nodes,1);
    
end

%% Get Mesh Displacements
disps=vecnorm(geom_temp.vertices-geom_temp.vertices_orig,2,2);
min_disps=quantile(disps,.05);
max_disps=quantile(disps,.95);
mean_disps=mean(disps);
median_disps=median(disps);
disps_ratio_table=table(min_disps,median_disps,max_disps,mean_disps)


%% get Dihedral Face Angle
edge_angles=getAllEdgeAngles(geom_temp.faces,geom_temp.vertices);
min_edge_angles=quantile(edge_angles,.05);
max_edge_angles=quantile(edge_angles,.95);
mean_edge_angles=mean(edge_angles);
median_edge_angles=median(edge_angles);
edge_angles_ratio_table=table(min_edge_angles,median_edge_angles,max_edge_angles,mean_edge_angles)
%% report data

min_aspect=quantile(aspects,.05);
max_aspect=quantile(aspects,.95);
mean_aspect=mean(aspects);
median_aspect=median(aspects);
aspect_ratio_table=table(min_aspect,median_aspect,max_aspect,mean_aspect)

min_skew=quantile(skewness,.05);
max_skew=quantile(skewness,.95);
mean_skew=mean(skewness);
median_skew=median(skewness);
skew_ratio_table=table(min_skew,median_skew,max_skew,mean_skew)

%% plot data
fig1=figure();
subplot(2,2,1);
% hist(aspects,100,'Normalization','probability')
cdfplot(aspects);
xlabel('Aspect Ratio (1 = ideal, INF = bad)')
title('Aspect Ratio');
% ytix = get(gca, 'YTick');
% num_vals=size(geom_temp.faces,1);
% set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/num_vals*100));


% fig2=figure();
subplot(2,2,2);
% hist(skewness,100,'Normalization','probability')
cdfplot(skewness);
xlabel('Skewness (0 = ideal, 1 = bad)')
title('Skewness');
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/num_vals*100))

%% average displacement
% fig3=figure()
subplot(2,2,3);
% hist(disps,100,'Normalization','probability')
cdfplot(disps)
xlabel('Displacement (mm)')
title('Displacement');
% ytix = get(gca, 'YTick');
% num_vals=size(geom_temp.vertices,1);
% set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/num_vals*100))

%% average angles
% fig4=figure()
subplot(2,2,4);
% hist(disps,100,'Normalization','probability')
cdfplot(edge_angles)
xlabel('Dihedral Face Angle (deg)')
title('Dihedral Face Angle');
% ytix = get(gca, 'YTick');
% num_vals=size(geom_temp.vertices,1);
% set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/num_vals*100))
