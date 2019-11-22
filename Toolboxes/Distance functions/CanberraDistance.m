function err = CanberraDistance(resp, targ),

err = sum(abs(resp(:) - targ(:))./(abs(resp(:)) + abs(targ(:))));