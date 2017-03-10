function plot_me( mytitle, cpu_nmb, cpu_times_relative, gpu_nmb, gpu_times_relative )
    hold on

    title(mytitle)
    
    set(gca,'fontsize',12);

    cpu_plot = plot(cpu_nmb, cpu_times_relative, 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_times_relative(1,:), 'bo', 'LineWidth', 2.0);

    gpu_plot = plot(gpu_nmb, gpu_times_relative, 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_times_relative(1,:), 'ro', 'LineWidth', 2.0);

%    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
%    ylabel('one iteration computing time [s]', 'Interpreter', 'latex', 'FontSize', 12)

%    h = legend([cpu_plot gpu_plot],'CPU','GPU');
%    set(h,'Interpreter','latex');
%    set(h, 'FontSize', 16);

    axis([1 max([cpu_nmb, gpu_nmb]) ...
      min([cpu_times_relative(1,:), gpu_times_relative(1,:)]) max([cpu_times_relative(1,:), gpu_times_relative(1,:)])])

    set(gca,'YScale','log');
    set(gca,'XScale','log');

    hold off

end

