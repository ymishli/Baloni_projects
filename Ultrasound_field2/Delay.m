function N_elements_delay = Delay(N_elements,c,pitch,focus_point, plot_for_debug, time_taps_generated_by_xdc_focus)
%This function calculates the delay for each element to get focus
%on focus_point
%   inputs:
%       - N_elements  - number of transducer elements
%       - c           - speed of sound
%       - pitch       - distance between two successive elements
%       - focus_point - the required point to be focus on
%   output:
%       - N_elements_delay - delay for each element 
if mod(N_elements,2) == 0 % even numbers of elements    
    positive_side = pitch*((0:N_elements/2-1)+1/2);
    negative_side = -flip(positive_side);
    elements_locations = [negative_side positive_side];
else % odd numbers of elements
    positive_side = pitch*((0:(N_elements-1)/2-1));
    negative_side = -positive_side;
    elements_locations = [negative_side 0 positive_side];
end
focus_point_norm = norm(focus_point);
element_focus_point_distance = sqrt((elements_locations - focus_point(1)).^2 + focus_point(3)^2);
N_elements_delay_offset = (focus_point_norm - element_focus_point_distance)/c;
N_elements_delay = N_elements_delay_offset;
if plot_for_debug
    figure('Name','Q2 - Section 3 - Eelements delay from geometrical calculation');
    stem(elements_locations,N_elements_delay,'filled', 'LineStyle','none');
    hold on;
    stem(elements_locations,time_taps_generated_by_xdc_focus','filled', 'LineStyle','none');
    hold off;
    legend('Created by Delay function','Created by xdc\_focus');
    title('Q2 - Section 3 - Elements delay from geometrical calculation');    
    xlabel('Elements Locations [mm]');ylabel('Elements Delay [sec]');
end
