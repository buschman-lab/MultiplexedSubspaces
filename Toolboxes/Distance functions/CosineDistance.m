function err = CosineDistance(resp, targ),

err = 1 - sum(resp(:).*targ(:))./(sqrt(sum(resp(:).^2))*sqrt(sum(targ(:).^2)));