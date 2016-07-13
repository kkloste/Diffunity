% FIEDLER CUT
%
% Display Fiedler cut on usps, and on netscience
%

clc;
addpath ../util; % for set_figure_size
load ../data/xy_1.mat ;

A = A|A';
A = A-diag(diag(A));
n = size(A,1);
d = full(sum(A,1))';
D = spdiags(d,0,n,n);
Dsqinv = spdiags(d.^(-1/2),0,n,n); 
P = colnormout(A);
L = D - A;
nL = Dsqinv*(L*Dsqinv);

[ V, lam] = eig(nL);
lam = diag(lam);

min_lam = min( lam( find(lam) ) );

fiedler = V(:,2);

fprintf( ' min_lam = %f ,   f^T*nL*f = %f \n' , min_lam, fiedler'*nL*fiedler );

[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,fiedler,'halfvol',true) ;





%% PLOT
cmappart = hot(9);
cmap = @() colormap(flipud(cmappart(2:6,:)));
[px,py] = gplot(A,xy);
gdraw = @() plot(px,py,'k-','LineWidth',0.4,'Color',0.55*[1,1,1]); 
sdraw = @() plot(xy(S,1),xy(S,2),'ko','MarkerSize',7);
bigmarkers = 15;
smallmarkers = 2;
figsize = [2.75,2.75];
figarea = [0,220,-600,-120];

epsmin = 1e-2;
alpha = 0.9;
vert = 211; % chosen as an illustrative example

%%
clf;

    f = x > 0;
    cla;
    gdraw(); hold on; 
    plot(xy((vert),1),xy((vert),2),'ko','MarkerSize',7);
    scatter(xy(f,1),xy(f,2),bigmarkers,log10(x(f)),'filled'); cmap();
    caxis([-3,0]);
    set_figure_size(figsize);
    axis(figarea);
 
    print( './images/fiedler_cut','-dpng','-r600','-painters');
