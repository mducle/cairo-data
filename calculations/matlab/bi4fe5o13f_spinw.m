function [bfof, spec] = bi4fe5o13f_spinw(J, phase, calc)
Jc1 = J(1);
Jc2 = J(2);
Jab1 = J(3);
Jab2 = J(4);
Jd = J(5);
%ff=11.6*2.5; Jc1 = 34/ff; Jc2 = 20/ff; Jab1 = 45/ff; Jab2 = 74/ff; Jd = 191/ff;

if nargin<2 || phase == 1
    S3 = 0.01;
else
    S3 = 2.5;
end

if nargin<3
    calc = 0;
end

if nargout > 1
    calc = 1
end

bfof = spinw();
bfof.genlattice('lat_const',[8.29950 8.29950 18.05730],'angled',[90 90 90],'sym','P 42/m b c');
bfof.addatom('r',[0.5  0.   0.0800],'S',2.5,'color',[0 0 255],'label','Fe3+');
bfof.addatom('r',[0.8515 0.8388 0], 'S',2.5,'color',[255 0 0],'label','Fe3+');
bfof.addatom('r',[0.5  0.   0.25  ],'S',S3,'color',[0 0 128],'label','Fe3+');

bfof.gencoupling('maxDistance',10);
bfof.addmatrix('value',Jc1,'label','Jc1','color','r');         bfof.addcoupling('mat','Jc1','bond',1);
bfof.addmatrix('value',Jc2,'label','Jc2','color',[128 0 0]);   bfof.addcoupling('mat','Jc2','bond',2);
bfof.addmatrix('value',Jab1,'label','Jab1','color','b');       bfof.addcoupling('mat','Jab1','bond',3);
bfof.addmatrix('value',Jab2,'label','Jab2','color',[0 255 0]); bfof.addcoupling('mat','Jab2','bond',4);
bfof.addmatrix('value',Jd,'label','Jd','color','k');           bfof.addcoupling('mat','Jd','bond',5);
if numel(J)==6
    bfof.addmatrix('value',diag([0 0 J(6)]),'label','D');      bfof.addaniso('D');
elseif numel(J)==7
    bfof.addmatrix('value',diag([0 0 J(6)]),'label','D1');     bfof.addaniso('D1', [1 2]);
    bfof.addmatrix('value',diag([0 0 J(7)]),'label','D2');     bfof.addaniso('D2', [3]);
elseif numel(J)==8
    bfof.addmatrix('value',diag([J(8) 0 J(6)]),'label','D1');  bfof.addaniso('D1', [1 2]);
    bfof.addmatrix('value',diag([J(8) 0 J(7)]),'label','D2');  bfof.addaniso('D2', [3]);
else
    bfof.addmatrix('value',diag([0 0 0.2]),'label','D');       bfof.addaniso('D');
end

S2a = -[4.05 -0.35 0];
S2b = -[-0.35 4.05 0];
S1a = [2.18 -2.53 0];
S1b = [2.53 2.18 0];
S = [S1b; S1a; S1a; S1b; S1b; S1a; S1a; S1b; S2a; -S2a; -S2b; S2b; S2b; -S2b; S2a; -S2a; -S1b; -S1a; -S1b; -S1a];
Sv=[S; -S; -S; S];

bfof.genmagstr('mode','direct','S',Sv','nExt',[2 2 1]);
plot(bfof, 'range', [0 0 0.4; 2 2 0.6]');

out = optmagsteep(bfof,'nRun',1000);
bfof = out.obj;

if calc
    q0 = [3 3 0];
    spec = bfof.spinwave({[-0.5 0 -0.5]+q0 [-0.5 0 0]+q0 [0 0 0]+q0 [0 -0.5 0]+q0 [-0.5 -0.5 0]+q0 [0 0 0]+q0 [-0.5 -0.5 0.5]+q0 [-0.5 -0.5 0]+q0 100},'hermit',false,'optmem',20);
    spec = sw_neutron(spec);
    try
        spech = sw_egrid(spec,'Evect',0:0.2:100.0); figure; sw_plotspec(spech,'mode','color','dE',1)
        title(sprintf('J_{c1}=%4.2f J_{c2}=%4.2f J_{ab1}=%4.2f J_{ab2}=%4.2f J_{d}=%4.2f', J));
    end
end
