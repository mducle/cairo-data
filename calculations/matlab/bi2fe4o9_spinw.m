function [bfo, spec] = bi2fe4o9_spinw(Jvals)

bfo = spinw;

lat_const = [7.9711 8.4391 6.0003];
angled = [90 90 90];

bfo.genlattice('lat_const', lat_const, 'angled', angled, 'spgr', 'P b a m');
bfo.addatom('r', [0 0 0.2583], 'S', 2.5, 'label', 'Fe3+', 'color', 'brown')
bfo.addatom('r', [0.1466 0.3361 0.5], 'S', 2.5, 'label', 'Fe3+', 'color', 'green')

bfo.gencoupling('maxDistance', 12);
%bondtable = bfo.table('bond', 1:9)

if nargin == 0
%Jvals = [2.3 2.0 3.6 5.1 13.5];
Jvals = [1.3 3.0 2.6 3.6 13.5];
end

for nn = 1:numel(Jvals)
    bfo.addmatrix('label', ['J' num2str(nn)], 'value', Jvals(nn));
    bfo.addcoupling('mat', ['J' num2str(nn)], 'bond', nn);
end

bfo.addmatrix('label', 'K', 'value', diag([0 0 0.021]))
bfo.addaniso('K');

S1 = -[2.1 1.3 0]';
S2 = [-0.9 2.3 0]';
S3 = [1.3 2.1 0]';
S4 = [2.3 -0.9 0]';
SS = [S1 S2 S1 S2 S3 S3 S4 S4];

bfo.genmagstr('mode', 'direct', 'nExt', [2 2 2], 'S', [SS -SS -SS SS -[SS -SS -SS SS]]);
bfo.energy

bfo.optmagsteep('nRun',50000)
bfo.energy
%bfo.magstr.S

plot(bfo, 'range', [2 2 2]);
%plot(bfo, 'range', [0.2 1.7; 0.2 1.7; 0 1]);

q0 = [3 3 0];
spec = bfo.spinwave({[-0.5 0 -0.5]+q0 [-0.5 0 0]+q0 [0 0 0]+q0 [0 -0.5 0]+q0 [-0.5 -0.5 0]+q0 [0 0 0]+q0 [-0.5 -0.5 0.5]+q0 [-0.5 -0.5 0]+q0 100},'hermit',false,'optmem',20);
%figure; sw_plotspec(spec,'mode','disp','imag',true,'colormap',[0 0 0],'colorbar',false);
spec = sw_neutron(spec);
spech = sw_egrid(spec,'Evect',0:0.2:100.0); figure; sw_plotspec(spech,'mode','color','dE',1)
%specl = sw_egrid(spec,'Evect',0:0.02:10.0); figure; sw_plotspec(specl,'mode','color','dE',1)
title(sprintf('J_{44}=%4.2f J_c=%4.2f J_{43}=%4.2f J_{43}''=%4.2f J_{33}=%4.2f', Jvals));
