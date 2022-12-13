function [Node] = GenerateNodes_omar(H,NH,V,NV)
%Omar's Code

Node=zeros((NH+1)*(NV+1),1);
if (NH==0) %Deals with zero horizontal elements
    Horz_Elem=0;
else
    Horz_Elem=H/NH; %Length of Horizontal Elements
end
    
 if (NV==0) %Deals with zero vertical elements
     Vert_Elem=0;
else
    Vert_Elem=V/NV; %Length of Vertical Elements
 end

        for i=1:NV+1 % i is for vertical
            for j=1:NH+1 % j is for horizontal
              Node((i-1)*(NH+1)+j,1)=(i-1)*(NH+1)+j;
              Node((i-1)*(NH+1)+j,2)=(j-1)*Horz_Elem;
              Node((i-1)*(NH+1)+j,3)=(i-1)*Vert_Elem;
            end
        end
end
                 