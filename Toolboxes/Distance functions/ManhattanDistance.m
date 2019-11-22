function err = ManhattanDistance(resp, targ),

err = sum(abs(resp(:) - targ(:)));