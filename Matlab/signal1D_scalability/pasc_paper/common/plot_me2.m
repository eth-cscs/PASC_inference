function plot_me2( mytitle, cpu_nmb, cpu_times_relative, gpu_nmb, gpu_times_relative )
    hold on

    title(mytitle)
    
    set(gca,'fontsize',12);

    cpu_plot = plot(cpu_nmb, cpu_times_relative(1)./cpu_times_relative, 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_times_relative(1)./cpu_times_relative, 'bo', 'LineWidth', 2.0);

    gpu_plot = plot(gpu_nmb, gpu_times_relative(1)./gpu_times_relative, 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_times_relative(1)./gpu_times_relative, 'ro', 'LineWidth', 2.0);

%    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
%    ylabel('one iteration computing time [s]', 'Interpreter', 'latex', 'FontSize', 12)

%    h = legend([cpu_plot gpu_plot],'CPU','GPU');
%    set(h,'Interpreter','latex');
%    set(h, 'FontSize', 16);
    plot(gpu_nmb,gpu_nmb,'k--');

    axis([1 max([cpu_nmb, gpu_nmb]) ...
      min([cpu_times_relative(1)./cpu_times_relative, gpu_times_relative(1)./gpu_times_relative]) max([cpu_times_relative(1)./cpu_times_relative, gpu_times_relative(1)./gpu_times_relative])])

%    set(gca,'YScale','log');
%    set(gca,'XScale','log');

    hold off

end

