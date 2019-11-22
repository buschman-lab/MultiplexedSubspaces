function err = MeanAbsError(resp, targ),

err = nanmean(abs(resp(:) - targ(:)));