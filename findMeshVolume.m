function [volume]=findMeshVolume(faces,vertices)
volume=0;
    for count_face=1:size(faces,1)
        face_vert=vertices(faces(count_face,:),:);
        x=face_vert(:,1);
        y=face_vert(:,2);
        z=face_vert(:,3);
        temp_vol=(1/6)*(-x(3)*y(2)*z(1)+x(2)*y(3)*z(1)+x(3)*y(1)*z(2)...
            -x(1)*y(3)*z(2)-x(2)*y(1)*z(3)+x(1)*y(2)*z(3));
        volume=volume+temp_vol;
    end

end