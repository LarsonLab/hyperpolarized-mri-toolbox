function p = rectangle_shape( x, y, z, Lx, Ly, Lz, x0, y0, z0 )
    p = (abs(x-x0)<Lx/2*ones(size(x))) .* (abs(y-y0)<Ly/2*ones(size(y))) .* (abs(z-z0)<Lz/2*ones(size(z)));
end
