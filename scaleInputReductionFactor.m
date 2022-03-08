function new_value=scaleInputReductionFactor(previous_value,scale_percent_factor)
    new_value=previous_value*scale_percent_factor;
    if new_value>=1
        new_value=1;
    end
end