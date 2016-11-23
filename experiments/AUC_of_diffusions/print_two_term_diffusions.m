function print_two_term_diffusions( which_comm, vert, varargin )
% print_two_term_diffusions( varargin )
%i
% p.addOptional('fname','senate');

p = inputParser;
p.addOptional('fname','senate');
p.parse(varargin{:});

fname = p.Results.fname;

save_dir = './results/';
image_dir = './images/';
addpath ../../util;

load( [save_dir, fname, '-', num2str(which_comm), '-', num2str(vert), '-two_term_diff.mat'] );

  stat_string = { 'cond', 'AUC' };
  Xs = diffusion_coeffs(:,1);
  linstyles = { '-', ':' };
  colors = [ 1, 0.5, 0.5; 0, 0, 1 ];
  for which_stat =1:2,
    clf
    for which_diff = 1:2,
      for deg_scale = 1:2,
          Ys = squeeze(diffusion_stats(which_diff, :, deg_scale, which_stat))';
        	plot( Xs,Ys, 'Color', colors(which_diff,:), 'LineStyle', linstyles{deg_scale} );
          %plot( Xs,Ys );
        	hold all;
      end
    end
    ylim([0,1]);
    xlim([0,max(Xs)]);

    title( sprintf( 'Graph: %s, community: %d , vert: %d', fname, which_comm, vert ) );
    xlabel('coeff 1');
    ylabel( stat_string{which_stat} );
    legend('full', 'full-n', 'push','push-n', 'location','Southeast');
    set_figure_size( [3,3] );
    print(gcf,[ image_dir, 'gen-diff-AUC-', stat_string{which_stat}, fname, '-', num2str(which_comm), '-', num2str(vert), '.png'],'-dpng');
  end



fprintf('Done printing  %s \n', fname);
