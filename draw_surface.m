function draw_surface(P, color)    

    x = P(1, :)';
    
    y = P(2, :)';
    
    z = P(3, :)';
    
%     [xi, yi] = meshgrid(linspace(min(x), max(x), 20), linspace(min(y), max(y), 20));
%     F = scatteredInterpolant(x,y,z);
%     zi = F(xi,yi);
    
%     surf(reshape(x, 20, 20),reshape(y, 20, 20),reshape(z, 20, 20), 'FaceAlpha',0.5)
%     surf(xi, yi, zi, 'FaceAlpha',0.5)
%     hold on
    plot3(x, y, z, ['o' color])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
