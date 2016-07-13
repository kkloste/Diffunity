% Fiedler plotting

load ./results/fiedler_usps.mat ;
load ../../data/usps3nn.mat ;

clf;

markersize = 2;

spy(A,'k', markersize);
title('sparsity pattern, USPS');
print( './images/usps3nn_spy','-dpng','-r600','-painters');

hold all

for j=1:length(bestset),
	vert = bestset(j);
	neighb = find(A(:,vert));
	neighb_in = intersect(neighb,bestset);
	scatter( neighb_in, vert*ones(length(neighb_in),1 ), markersize, 'b', 'filled' );
	scatter( vert*ones(length(neighb_in),1 ), neighb_in, markersize, 'b', 'filled' );

	neighb_out = setdiff(neighb, neighb_in);
	scatter( neighb_out, vert*ones(length(neighb_out), 1 ), markersize, 'r', 'filled' );
	scatter( vert*ones(length(neighb_out), 1 ), neighb_out, markersize, 'r', 'filled' );

end

title('fiedler cut, USPS')
print( './images/usps3nn_fiedler','-dpng','-r600','-painters');


