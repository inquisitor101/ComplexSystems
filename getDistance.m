function dist = getDistance(i, Cx, Cy, N)

dist = zeros(N-1, 1); % initialize distance

for j=1:i-1
    x = Cx(i) - Cx(j);
    y = Cy(i) - Cy(j);
    dist(j) = sqrt( x*x + y*y );
end

for j=i+1:N
    x = Cx(i) - Cx(j);
    y = Cy(i) - Cy(j);
    dist(j) = sqrt( x*x + y*y );
end