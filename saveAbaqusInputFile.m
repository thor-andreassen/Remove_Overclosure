%% clearing
clear
close all
clc

%% load data
load('muscle_geom_orig.mat');

%% create filename
filename='S192803L_Femur_DEFORM.inp';

%% data
fid=fopen(filename,'w+');
fprintf(fid,'*NODE,NSET=FEMUR_NODES\r\n');
geom=geom2;

for countnode=1:size(geom.new_nodes,1)
    fprintf(fid,'%d, %.6f, %.6f, %.6f\r\n',geom.new_nodes(countnode,1),...
        geom.new_nodes(countnode,2),geom.new_nodes(countnode,3),geom.new_nodes(countnode,4));
end

fprintf(fid,'**\r\n**\r\n**\r\n');
if size(geom.original_elements,2)==4
    elem_type='R3D3';
elseif size(geom.original_elements,2)==5
    elem_type='C3D4';
elseif size(geom.original_elements,2)==6
    elem_type='C3D5';
elseif size(geom.original_elements,2)==9
    elem_type='C3D8R';
end
fprintf(fid,['*ELEMENT, TYPE=',elem_type,', ELSET=FEMUR_ELEMS\r\n']);
for countelems=1:size(geom.original_elements,1)
    fprintf(fid,'%d, ',geom.original_elements(countelems,1));
    for count_node=2:size(geom.original_elements,2)
        if count_node==size(geom.original_elements,2)
            fprintf(fid,'%d\r\n',geom.original_elements(countelems,count_node));
        else
            fprintf(fid,'%d, ',geom.original_elements(countelems,count_node));
        end
    end
end
fprintf(fid,'**\r\n**\r\n**\r\n');
fclose(fid);