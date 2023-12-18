function derivatives = calculate_derivative(~, states, m1, m2, k, c, l, g)
    x1 = states(1);
    y1 = states(2);
    vx1 = states(3);
    vy1 = states(4);
    x2 = states(5);
    y2 = states(6);
    vx2 = states(7);
    vy2 = states(8);
    
    distance = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    spring_force = k * (distance - l);
    damper_force_x = c * (vx2 - vx1);
    damper_force_y = c * (vy2 - vy1);
    
    ax1 = (spring_force * (x2 - x1)) / (m1 * distance) + damper_force_x / m1;
    ay1 = (spring_force * (y2 - y1)) / (m1 * distance) + damper_force_y / m1 - g;
    
    ax2 = (spring_force * (x1 - x2)) / (m2 * distance) - damper_force_x / m2;
    ay2 = (spring_force * (y1 - y2)) / (m2 * distance) - damper_force_y / m2 - g;
    
    derivatives = [vx1; vy1; ax1; ay1; vx2; vy2; ax2; ay2];
end



