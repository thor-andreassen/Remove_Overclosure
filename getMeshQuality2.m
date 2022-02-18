function [skewness, aspect_ratio]=getMeshQuality2(nodes,is2D)
    if is2D
        if size(nodes,1)==3
            tri_face=1;
            edges=[1,2;2,3;3,1];
            theta_ideal=60;
        else
            tri_face=0;
            edges=[1,2;2,3;3,4;4,1];
            theta_ideal=90;
        end
    else
        if size(nodes,1)==4
            tri_face=1;
            edges=[1,2;2,3;3,1;1,4;2,4;3,4;];
        else
            tri_face=0;
            edges=[1,2;2,3;3,4;4,1;5,6;6,7;7,8;8,5;1,5;2,6;3,7;4,8];
        end
    end
    
    edge_lengths=zeros(size(edges,1),1);
    angles=edge_lengths;
    current_edge_node=zeros(2,3);
    edge_vectors=zeros(size(edges,1),3);
    
    for count_edge=1:size(edges,1)
       current_edge_node(1,:)=nodes(edges(count_edge,1),:);
       current_edge_node(2,:)=nodes(edges(count_edge,2),:);
       edge_lengths(count_edge)=norm(current_edge_node(2,:)-current_edge_node(1,:));
       edge_vectors(count_edge,:)=(current_edge_node(2,:)-current_edge_node(1,:))/...
       norm(current_edge_node(2,:)-current_edge_node(1,:));
    end
    
    for count_edge=1:size(edges,1)
       vec1=-edge_vectors(edges(count_edge,1),:);
       vec2=edge_vectors(edges(count_edge,2),:);
       angles(count_edge)=acos(dot(vec1,vec2));
       angles(count_edge)=abs(angles(count_edge)*(180/pi));
    end
    
    skewness=max([(max(angles)-theta_ideal)/(180-theta_ideal),...
        (theta_ideal-min(angles))/(theta_ideal)]);
    
    aspect_ratio=max(edge_lengths)/min(edge_lengths);
end