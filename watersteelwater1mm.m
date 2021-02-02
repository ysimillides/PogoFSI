clear

%Ly is the layer height (1mm+ 10mm either side)
Ly = 0.021 %0.006
Lx = 0.6 %0.025

dx = 5e-5;
dy = 5e-5;

nx = round((Lx/dx)+1)
ny = round((Ly/dy)+1)

cx = ((dx/2)*nx)-(dx/2);
cy = ((dy/2)*ny)-(dy/2);

m = genGrid2D(nx,ny,dx,dy,cx,cy);                                       
%imports mode shapes from Disperse txt file, Uxreal/Uzreal are the ones
%we're interested in
[Position,Uxreal,Uximag,Uzreal,Uzimag,Uyreal,Uyimag] = importfile4('D:\OneDrive - Imperial College London\watersteelwater.txt',2, 52)

%%
nNodes = size(m.nodePos,2);
nEls = size(m.elNodes,2);
%% 

%this is steel
E = 216.90e9;
nu = 0.2865;
G = E/(2*(1+nu))
rho = 7932;

%this is brass
%Eb = 108.41e9;
%nub = 1/3;
%Gb = Eb/(2*(1+nub));
%rhob = 8400

%this is water
rhow = 1000
% Ew= rhow*1238*1238
cw = 1500;
nuw = 0%0.5-aa

%this is aluminum.
% Ea = 70.758e9;
% nua = 0.3375;
% Ga= E/(2*(1+nu))
% rhoa = 2700;

%this is magnesium
% Em =41.310e9;
% num = 0.3061;
% Gm = E/(2*(1+nu))
% rhom = 1700;


%more conservative to be consistent with fluid
courant = 0.3;

freq =1000e3;
nCyc =9;
cL = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu)));
cSh = sqrt(E/(2*rho*(1+nu)));

endTime = Lx*2/cL;                          % Model run time 
CFL = courant                                
dt = CFL*dx/cL                         % Stable Time with Courant number CFL

lambda = cL/freq;
lambdaSh = cSh/freq;

cMax = cL;

nt = round(endTime / dt) 
m.nt = nt;
m.dt = dt;

nBound = 60;
alpha = 0;
m.matTypeRefs = ones(nEls,1);
px = m.nodePos(1,:);
py = m.nodePos(2,:);
%%
[ex, ey] = getElCents(m);

%define the "layer of interest"
bottom = (Ly-0.001)/2;
top = bottom+0.001;
middle = (top+bottom)/2

%find elements either side of layer of interest
waterEls = find(ey <bottom | ey>top);

m.elTypes = cell(2,1);
m.matTypeRefs(waterEls) = 2 
m.elTypeRefs(waterEls) = 2
m.matTypes = cell(2,1);

% m.elTypeRefs - which of the element types each element refers to, length nEls
m.elTypes{1}.name = 'CPE4R'
m.elTypes{2}.name = 'CPE4R'


m.elTypes{1}.paramsType = 0;
m.elTypes{2}.paramsType = 0;

m.matTypes{1}.paramsType = 0;
m.matTypes{1}.paramValues = [E, nu, rho];


%define fluid parameters
m.matTypes{2}.paramsType = 5;
m.matTypes{2}.paramValues = [rhow,cw,nuw];
%% 
%need to generate two signals here, (use sin/cos to make them out of phase
%by 90 degrees)
m = genPogoHannSignal(m,nCyc,freq,1,1,0);
m = genPogoHannSignal(m,nCyc,freq,2,1,1);

siglength = (1/freq)*nCyc


%%
measnode = []
mylocation = []

nsrc = 21%/40 %changed to 40 overnight

%sigs{1} refers to first dof
%sigs{2} refers to second dof

m.shots{1}.sigs{1}.nodeSpec = zeros(nsrc,1);
m.shots{1}.sigs{1}.dofSpec = zeros(nsrc,1);
m.shots{1}.sigs{2}.nodeSpec = zeros(nsrc,1);
m.shots{1}.sigs{2}.dofSpec = zeros(nsrc,1);

%find locations of sources from disperse data
ylocs = linspace(bottom,top,nsrc);
sourcelocs=[];
sourcepos=linspace(1,51,nsrc);
sourcepos = round(sourcepos);

% find node locations of sources
for i=1:nsrc
   mylocation = find(abs(py-ylocs(i))<0.5*dy &abs(px-0-0.0)<0.1*dx);
   %mylocation = find(abs(py-0.003-i*spacing*dy)<0.5*dy&abs(px-0.04)<0.5*dx);
   sourcelocs = [sourcelocs mylocation];
   m.shots{1}.sigs{1}.nodeSpec(i) = mylocation;
   m.shots{1}.sigs{2}.nodeSpec(i) = mylocation ;
   m.shots{1}.sigs{1}.dofSpec(i) = 1;
   m.shots{1}.sigs{2}.dofSpec(i) = 2;

end

%set Amplitudes based on Disperse data
sigAmpsx = Uxreal(sourcepos);
sigAmpsy = Uzreal(sourcepos);
normc = max(abs(sigAmpsy))
normd = max(abs(sigAmpsx))

%normalize this
sigAmpsx = sigAmpsx / normc;
sigAmpsy = sigAmpsy / normc;

m.shots{1}.sigs{1}.sigAmps = sigAmpsx;
m.shots{1}.sigs{1}.sigType = 0;

m.shots{1}.sigs{2}.sigAmps = sigAmpsy;
m.shots{1}.sigs{2}.sigType = 0;

m.measFreq = 20;
m.measStart = 1;

%measuring locations, along midpoint and top boundary
m.measSets{1}.name = 'midpointx';
m.measSets{2}.name = 'midpointy';
m.measSets{3}.name = 'topx';
m.measSets{4}.name = 'topy';


measnode = [];
measnodetop = [];
spacing = 4;

%find node locations of the measuring locations
for i=1:10000
   mylocation = find(abs(py-middle)<0.1*dy &abs(px-0.-i*spacing*dx)<0.1*dx);
   mylocation2 = find(abs(py-top)<0.1*dy &abs(px-0.-i*spacing*dx)<0.1*dx);
   %mylocation = find(abs(py-0.003-i*spacing*dy)<0.5*dy&abs(px-0.04)<0.5*dx);
   measnode = [measnode mylocation];
   measnodetop = [measnodetop mylocation2];
end
m.measSets{1}.measNodes = measnode
m.measSets{2}.measNodes = measnode
m.measSets{3}.measNodes = measnodetop
m.measSets{4}.measNodes = measnodetop

m.measSets{1}.measDof = (ones(1,length(measnode)));
m.measSets{2}.measDof = (ones(1,length(measnode))*2);
m.measSets{3}.measDof = (ones(1,length(measnodetop)));
m.measSets{4}.measDof = (ones(1,length(measnodetop))*2);

nt = m.nt
gap = round(nt/60);
if gap < 1
    gap = 1;
end

%not necessary at present
m.fieldStoreIncs = 1:gap:nt;

% %% Saving output
savePogoInp('watersteelwater1mm.pogo-inp',m);
%savePogoInp('C:\Users\yiann\OneDrive - Imperial College London\results\fluidlambwave9MARCH.pogo-inp',m);
