
%% test
clear CData CData_vec
face_list1=getHexorTetFaces(geom1.faces);
face_list2=getHexorTetFaces(geom2.faces);
CData=geom1.vertices-geom1_o.vertices_orig;
for counti=1:size(CData,1)
      CData_vec(counti)=norm(CData(counti,:));
end
CData_vec(counti)=norm(CData(counti,:));
% [~,CData_vec]=sort(CData_vec);
% color_set=jet(size(CData,1));
% color_set=color_set(CData_vec,:);
figure()
% p1=patch('Faces',face_list1,'Vertices',geom1.vertices,'FaceColor','interp','EdgeAlpha',.1,'FaceVertexCData',color_set);
for counti=1:size(face_list1,1)
        current_face=face_list1(counti,:);
        patch(geom1.vertices(current_face,1),geom1.vertices(current_face,2),geom1.vertices(current_face,3),CData_vec(current_face));
        hold on
end
clear CData CData_vec
CData=geom2.vertices-geom2_o.vertices_orig;
for counti=1:size(CData,1)
      CData_vec(counti)=norm(CData(counti,:));
end
CData_vec(counti)=norm(CData(counti,:));
% [~,CData_vec]=sort(CData_vec);
% color_set=jet(size(CData,1));
% color_set=color_set(CData_vec,:);
% p2=patch('Faces',face_list2,'Vertices',geom2.vertices,'FaceColor','interp','EdgeAlpha',.1,'FaceVertexCData',color_set);
for counti=1:size(face_list2,1)
      current_face=face_list2(counti,:);
      patch(geom2.vertices(current_face,1),geom2.vertices(current_face,2),geom2.vertices(current_face,3),CData_vec(current_face));
      hold on
end
colormap jet
colorbar
