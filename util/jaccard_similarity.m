function [sim, intersz, unionsz] = jacccard_similarity( set1, set2 )
% [jaccard_sim, intersection_size, union_size] = jacccard_similarity( set1, set2 )

intersz = length( intersection(set1, set2 ) );
unionsz = length( union(set1, set2 ) );

sim = 1.0;

if unionsz == 0, 
	return;
end

sim = intersz/unionsz;
