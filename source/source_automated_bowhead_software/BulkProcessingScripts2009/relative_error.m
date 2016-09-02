
function err=relative_error(test_feature,anchor_feature)
    err=abs(test_feature-anchor_feature)./anchor_feature;
end

