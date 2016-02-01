clc; clear all; close all
n_particles = 20;
box_size = 4.0;
nside = 4;

x = linspace(-box_size/2,box_size/2,nside+1);
hold on
for k=1:numel(x)
    plot([x(k) x(k)], [x(1) x(end)])
    plot([x(1) x(end)],[x(k) x(k)])
end

number_of_boxes = nside^2;
for i=1:n_particles
    z_x(i) = box_size*rand-box_size/2.0;
    z_y(i) = box_size*rand-box_size/2.0;
end

for k=1:n_particles
    plot(z_x(k),z_y(k),'*');
end

nparticles_in_box = zeros(number_of_boxes,1);

for j = 1:n_particles
    box_x = floor(nside*(z_x(j)+box_size/2.0)/box_size);
    box_y = floor(nside*(z_y(j)+box_size/2.0)/box_size);
    if(box_x < 0) box_x = 0; end
    if(box_x >= nside) box_x = nside-1; end
    if(box_y < 0) box_y = 0; end
    if(box_y >= nside) box_y = nside-1; end
    in_box(j) =  box_y*nside + box_x;
    nparticles_in_box(in_box(j)+1) = nparticles_in_box(in_box(j)+1) +1;
end
nparticles_in_box = reshape(nparticles_in_box,nside,nside);
nparticles_in_box = rot90(nparticles_in_box)