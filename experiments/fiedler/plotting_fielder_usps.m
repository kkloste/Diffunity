% Fiedler plotting

load ./results/fiedler_usps.mat ;
load ../../data/usps3nn.mat ;

clf;

markersize = 6;

spy(A,'k', markersize);
title('sparsity pattern, USPS');
print( './images/usps3nn_spy','-dpng','-r600','-painters');

%for j=1:length(bestset),
for j=1:20,
	vert = bestset(j);
	neighb = A(:,vert);
	neighb_in = intersect(neighb,bestset);
	scatter( neighb_in, vert*ones(length(neighb_in),1 ), 'b', markersize, 'filled' );
	scatter( vert*ones(length(neighb_in),1 ), neighb_in, 'b', markersize, 'filled' );

	neighb_out = setdiff(neighb, neighb_in);
	scatter( neighb_out, vert*ones(length(neighb_out), 1 ), 'r', markersize, 'filled' );
	scatter( vert*ones(length(neighb_out), 1 ), neighb_out, 'r', markersize, 'filled' );

end

title('fiedler cut, USPS')
print( './images/usps3nn_fiedler','-dpng','-r600','-painters');


