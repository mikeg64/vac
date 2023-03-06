xirr=[0 3 4; 5 0 5]';
yirr=[5 4 3; 0 0 5]';
nx=40; ny=40;

triangulate(xirr,yirr,nx,ny);

zirr=[3. 3. 3.; 3. 2 1]';
zreg=interpolate(zirr);
clf reset; figure(gcf);
mesh(zreg);
