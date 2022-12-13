function [elem] = trail_elements(NH,NV)  
% Omar's code

for m=1:NV+1
    for n=1:NH
elem((m-1)*(NH)+n,1)=(m-1)*(NH)+n;
elem((m-1)*(NH)+n,2)=(m-1)*(NH)+n+m-1;
elem((m-1)*(NH)+n,3)=(m-1)*(NH)+n+m;
  if(m==NV+1 && n==NH)
      break
  end
    end
end
    
      for x=1:NH+1
          for y=1:NV
              elem((x-1)*(NV)+y+(NV+1)*NH,1)=(x-1)*(NV)+y+(NV+1)*NH;
              elem((x-1)*(NV)+y+(NV+1)*NH,2)=(x-1)*(NV)+y+(NV+1)*NH-(NV*NH+NH)-(NV-1)*(x-1)+(y-1)*(NH+1)-(y-1);
              elem((x-1)*(NV)+y+(NV+1)*NH,3)=(x-1)*(NV)+y+(NV+1)*NH-(NV*NH+NH)-(NV-1)*(x-1)+(y-1)*(NH+1)-(y-1)+NH+1;
          end
     end

end

