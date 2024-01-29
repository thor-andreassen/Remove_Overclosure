function writeAbaqusInput(file_name,nodes,elements,params)
    elem_type=params.elem_type;
    nset_name=params.nset_name;
    elset_name=params.elset_name;
    
    fid=fopen(file_name,'w+');
    
    fprintf(fid,'**\r\n');
    fprintf(fid,'**\r\n');
    fprintf(fid,'**\r\n');
    fprintf(fid,'*NODE, NSET=%s\r\n',nset_name);
    for count_node=1:size(nodes,1)
        fprintf(fid,'%d, %12f, %12f, %12f\r\n',...
            nodes(count_node,1),...
            nodes(count_node,2),...
            nodes(count_node,3),...
            nodes(count_node,4));
    
    end
    
    
    fprintf(fid,'**\r\n');
    fprintf(fid,'**\r\n');
    fprintf(fid,'*ELEMENT,TYPE=%s,ELSET=%s\r\n',...
        elem_type,nset_name);
    
    for count_row=1:size(elements,1)
        for count_col=1:size(elements,2)
            if count_col<size(elements,2)
                fprintf(fid,'%d, ',elements(count_row,count_col));
            else
                fprintf(fid,'%d\r\n',elements(count_row,count_col));
            end
        end
    end
fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');
fclose(fid);
end
