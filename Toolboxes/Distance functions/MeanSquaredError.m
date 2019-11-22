function err = MeanSquaredError(resp, target),

err = nanmean((resp(:) - target(:)).^2);