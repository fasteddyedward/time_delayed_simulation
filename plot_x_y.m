function plot_x_y(movie_name)
load([movie_name,'.mat'],'x','y','N','time');
for i=1:N
    figure(100+i)
    hold on
    plot(time,x(i,1:end-1))
    plot(time,y(i,1:end-1))
    title(['Position of particle ',num2str(i)])
    legend('x','y')
    xlabel('time (ms)')
    ylabel('position (mm)')
end
end