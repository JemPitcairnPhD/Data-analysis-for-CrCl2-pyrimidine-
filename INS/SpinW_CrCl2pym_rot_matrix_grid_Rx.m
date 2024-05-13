%% define nuclear lattice
CrPymCl2 = spinw;
CrPymCl2.genlattice('lat_const',[3.6688 12.103 7.0628],'angled',[90 94.236 90],'spgr','P 21/m')
CrPymCl2.addatom('r',[0 0 0],'S', 2,'label','Cr2','color','b')
CrPymCl2.table('atom')
plot(CrPymCl2)
swplot.zoom(1.5)


CrPymCl2.gencoupling('maxDistance',7.5)
CrPymCl2.table('bond',[])
J1 = 1.13 %AFM
J2 = -0.1 % FM
J3 = -0.01 % FM
D1 = -0.11 %Ising
CrPymCl2.addmatrix('label','J1','value',J1,'color','red')
CrPymCl2.addmatrix('label','J2','value',J2,'color','green')
CrPymCl2.addmatrix('label','J3','value',J3,'color','pink')
CrPymCl2.addcoupling('mat','J1','bond',1)
CrPymCl2.addcoupling('mat','J2','bond',2)
CrPymCl2.addcoupling('mat','J3','bond',3)

vec = [1.82 0 -0.75]
v = [-1 0 0]% starting axis vector - cannot be off axis, because then the matrix is incorrect
mag = dot(vec,vec)^0.5 %get magnitude for the end
vec = vec/mag %normalise final unit vector
c = cross(v,vec); % cross product
e = c/(dot(c,c)^0.5); %normalised cross 
d = dot(v,vec); %dot product
si = dot(c,c)^0.5; %sin of the angle
K = [0 -1*e(3) e(2); e(3) 0 -1*e(1); -1*e(2) e(1) 0]; %'cross product matrix'

rot = diag([1 1 1])+si*K+(1-d)*K*K; %axis angle formula
ani=rot*diag(v)*inv(rot);

d_mat=ani.*-1;
CrPymCl2.addmatrix('value',d_mat.*D1,'label','D','color','r');
CrPymCl2.addaniso('D')

plot(CrPymCl2,'range',[2 2 1.0])

% Define Magnetic lattice
CrPymCl2.genmagstr('mode','helical','k',[0.5 0.0 0.0],'n',[0 1 0 ], 'S',[1.8; 0; -0.7],'nExt',[2 1 1]);

disp('Magnetic structure:')
CrPymCl2.table('mag')

%CrPymCl2.energy
plot(CrPymCl2,'range',[2 2 1])


% Optimise magnetic lattice
CrPymCl2.optmagsteep()
%CrPymCl2.table('mag')
CrPymCl2.energy
%%
plot(CrPymCl2,'range',[2 2 1])
disp('Magnetic structure:')
CrPymCl2.table('mag')

% Calculate spin waves
maxE = 9.71;
minE = -1.21;
dE=0.03035;
Es = minE+dE/2:dE:maxE;

dQ=0.02112;
minQ = 0.01;
maxQ = 4.68;
Qs = minQ+dQ/2:dQ:maxQ;
handle = figure;
CrPymCl2_spec = CrPymCl2.powspec(Qs,'Evect',Es,'nRand',1e4,'T',1.5,'formfact',true);

%%
CrPymCl2_spec = sw_instrument(CrPymCl2_spec,'dE',0.3,'dQ',0.1);% specify broadening here
name = (['CrpymCl2_',num2str(J1),'_',num2str(J2),'_',num2str(J3),'_',num2str(D1),'_PhD_symposium.dat'])

writematrix(CrPymCl2_spec.swConv, name,'Delimiter','tab') 
handle = figure
name = (['grid_CrPymCl2_',num2str(J1),'_',num2str(J2),'_',num2str(J3),'_',num2str(D1),'_Rx.png'])
sw_plotspec(CrPymCl2_spec,'mode',3);
axis([0 4 0 10])
caxis([0 0.1])
colormap(cm_viridis)
colorbar('off');
%saveas(handle,name)