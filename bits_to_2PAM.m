function [array] = bits_to_2PAM(b)
    
    for i = 1:length(b)
        if b(i) == 0
            array(i) = 1;
        else
            array(i) = -1;
        end
    end
end