function [] = animateThis(maxInformed, N, L, t, h, Cx, Cy, Xc, Yc, pauseTime, isPeriodic)


for i=1:maxInformed
   idx_x = abs( Cx(i, t+1) - Cx(i, t) );
   idx_y = abs( Cy(i, t+1) - Cy(i, t) );
   if idx_x < 0.5*L && idx_y < 0.5*L
       plot([Cx(i, t), Cx(i, t+1)], ...
            [Cy(i, t), Cy(i, t+1)], ...
            'r-','markersize',4)
       if isPeriodic
           axis([0 L 0 L]);
       end
       hold on
       plot(Cx(i, t+1), Cy(i, t+1), 'r.', ...
               'markersize', 10)
       xlabel('X position')
       ylabel('Y position')
       title(['time: ',num2str(t), ...
       '   informed (red): ',num2str(maxInformed), ... 
       ' h: ',num2str(rad2deg(h(t+1)))]);
   end
end
for i=maxInformed+1:N
   idx_x = abs( Cx(i, t+1) - Cx(i, t) );
   idx_y = abs( Cy(i, t+1) - Cy(i, t) );

   if idx_x < 0.5*L && idx_y < 0.5*L
       plot([Cx(i, t), Cx(i, t+1)], ...
            [Cy(i, t), Cy(i, t+1)], ...
            'b-','markersize',4)
       if isPeriodic
           axis([0 L 0 L]);
       end
       hold on
       plot(Cx(i, t+1), Cy(i, t+1), 'b.', ...
            'markersize', 10)
       xlabel('X position')
       ylabel('Y position')
       title(['time: ',num2str(t), ...
       '   informed (red): ',num2str(maxInformed), ...
       ' h: ', num2str(rad2deg(h(t+1))) ]);
   end
end

quiver(Xc(t), Yc(t), cos(h(t+1)), sin(h(t+1)), 'k', 'MaxHeadSize', 1.0 )

% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% debugging here (remove) ------------------------------
legend(['direction: ', num2str(rad2deg(h(t+1)))]);    
% debugging here (remove) ------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
getframe();
pause(pauseTime); hold off

