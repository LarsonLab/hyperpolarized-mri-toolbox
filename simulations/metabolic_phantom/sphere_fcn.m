function p = sphere_fcn(x, y, z, radius, x0, y0, z0 )
    p = ((x-x0).^2 + (y-y0).^2 + (z-z0).^2 < (radius)^2); 
end