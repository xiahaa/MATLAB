paper track:

| [On-board real-time pose estimation for uavs using deformable visual contour registration](javascript:void(0))A Amor-Martinez, A Ruiz, F Moreno-Noguer, A Sanfeliu2014 IEEE International Conference on Robotics and Automation (ICRA), 2595-2601 | [19](https://scholar.google.com/scholar?oi=bibs&hl=zh-CN&cites=5718035092537380329) | 2014 |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ---- |
| [Planar PØP: Feature-less pose estimation with applications in UAV localization](javascript:void(0))A Amor-Martinez, A Santamaria-Navarro, F Herrero, A Ruiz, A Sanfeliu2016 IEEE International Symposium on Safety, Security, and Rescue Robotics … | [1](https://scholar.google.com/scholar?oi=bibs&hl=zh-CN&cites=8047553645920362255) | 2016 |
| [Precise Localization for Aerial Inspection Using Augmented Reality Markers](javascript:void(0))A Amor-Martinez, A Ruiz, F Moreno-Noguer, A SanfeliuAerial Robotic Manipulation, 249-259 |                                                              |      |



(7) (8) is actually the integration (sum) of all points on the segments of their jacobian relative to the motion parameters. for example

for rotation 2D

R=[cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];

a = 0

R=I

R(da)=[1 -da 0;da 1 0;0 0 1];

then the jacobian is like

x-day-x/da=-y

y+dax-y/da=x

so for rotation, it is M(-1 0 1 1 1 0)

for homography, it does like

x/(dhx+1)-x/dh=-x^2

y/(dhx+1)-y/dh=-xy

so M(-1,2,0,-1,1,1)



for rigid motion, this can be done with projection and se3 using chain rule.