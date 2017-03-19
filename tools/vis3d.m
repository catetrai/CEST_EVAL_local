function vis3d(x)
% visualize 3d array

figure;
slice(x, round(size(x,1)/2), round(size(x,2)/2), round(size(x,3)/2));
xlabel('x'); ylabel('y'); zlabel('z');
grid on, shading interp, colormap gray