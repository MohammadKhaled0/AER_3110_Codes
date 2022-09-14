'(ch10p5) Example 10.11' % Display label.
clf % Clear graph.
numg=50; % Define numerator of G(s).
deng=poly([0 -3 -6]); % Define denominator of G(s).
'G(s)' % Display label.
G=tf(numg,deng) % Create and display G(s).
'T(s)' % Display label.
T=feedback(G,1) % Find and display closed-loop
% transfer function.
bode(T) % Make a Bode plot.
grid on % Turn on the grid for the plots.
title('Closed-Loop Frequency Response')
% Add a title to the Bode plot.
pause
nyquist(T) % Make a Nyquist diagram.
title('Closed-Loop Frequency Response')
% Add a title to the Nyquist
% diagram.
pause
ch10p6