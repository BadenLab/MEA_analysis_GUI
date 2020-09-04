function out = get_fileformat (in) 
%Finds and returns the part of the string behind the last dot in the string

    point_position = strfind(in,'.');
    
    point_position = point_position(end);
    
    out = in(point_position:end);
    
    

end