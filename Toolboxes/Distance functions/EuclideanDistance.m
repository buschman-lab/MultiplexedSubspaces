function err = EuclideanDistance(x, y)
 err = sqrt(sum((x(:) - y(:)).^2));

