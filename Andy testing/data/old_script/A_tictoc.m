clear;
clc;
t = zeros(1,100);
for n = 1:100
    A = rand(n,n);
    b = rand(n,1);
    tic;
    x = A\b;
    t(n) = toc;
    disp( num2str( toc ) );
end
plot(t)