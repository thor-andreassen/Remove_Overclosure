clear
close all
clc


elem_type='C3D8R';
nset_name='FEMUR_NODES';
elset_name='FEMUR_ELEMS';





node_nums=[100:199]';
nodes=rand(100,3);
nodes=[node_nums,nodes];


elements=[1000:1999]';


elem_num=8;

elem_rand=zeros(length(elements),elem_num);
for count_row=1:length(elements)
    elem_rand(count_row,:)=randperm(length(node_nums),elem_num)...
        +min(node_nums);
end


elements=[elements,elem_rand];

params.elem_type='C3D8R';
params.nset_name='NODES_FEMUR';
params.elset_name='ELEMENTS_FEMUR';

writeAbaqusInput('test_input.inp',nodes,elements,params)
% fid=fopen('test_input.txt','w+');
% 
% fprintf(fid,'**\r\n');
% fprintf(fid,'**\r\n');
% fprintf(fid,'*NODE, NSET=%s\r\n',nset_name);
% for count_node=1:length(node_nums)
%     fprintf(fid,'%d, %12f, %12f, %12f\r\n',...
%         node_nums(count_node),...
%         nodes(count_node,1),...
%         nodes(count_node,2),...
%         nodes(count_node,3));
% 
% end
% 
% 
% fprintf(fid,'**\r\n');
% fprintf(fid,'**\r\n');
% fprintf(fid,'*ELEMENT, ELSET=%s, TYPE=%s\r\n',...
%     nset_name,elem_type);
% 
% for count_row=1:size(elements,1)
%     for count_col=1:size(elements,2)
%         if count_col<size(elements,2)
%             fprintf(fid,'%d, ',elements(count_row,count_col));
%         else
%             fprintf(fid,'%d\r\n ',elements(count_row,count_col));
%         end
%     end
% end
% 
% fclose(fid)
% 
% patch('Faces',elements(:,2:4)-100,'Vertices',nodes(:,2:end),...
%     'FaceColor','r','EdgeAlpha',.3)
% 
%    
   
   
    