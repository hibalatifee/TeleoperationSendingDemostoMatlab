function [ Data ] = FixPolarRange( Data,min_val,max_val )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% while(1)
%     [val1 ind1]=find(Data<min_val);
%     Data(ind1)=Data(ind1)+2*pi;
%     [val2 ind2]=find(Data>max_val);
%     Data(ind2)=Data(ind2)-2*pi;
%     if(isempty(ind1)&&isempty(ind2) )
%        break; 
%     end
% end    
while(1)
    Range_Val = max_val - min_val;
    Avg_Val = min_val + (max_val - min_val)/2;
    [val1 ind1]=find(Data<min_val);
    Data(ind1) = Data(ind1) + Range_Val * floor(abs(Data(ind1)-max_val)/Range_Val);
    [val2 ind2]=find(Data>max_val);
    Data(ind2) = Data(ind2) - Range_Val * floor(abs(Data(ind2)-min_val)/Range_Val);
     if(isempty(ind1)&&isempty(ind2) )
        break; 
     end
end
end

