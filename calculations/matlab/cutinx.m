function dat = cutinx(inx, ax, ll, ul)

if ax == 'x'
    id = find((inx.x>=ll) .* (inx.x<ul));
    yy = inx.z(:, id)';
    ee = inx.e(id, :);
    dat.x = inx.y;
else
    id = find((inx.y>=ll) .* (inx.y<ul));
    yy = inx.z(id, :);
    ee = inx.e(:, id)';
    dat.x = inx.x;
end

dat.y = dat.x * NaN;
dat.e = dat.x * NaN;

for ii = 1:size(yy, 2)
    y0 = yy(:, ii);
    e0 = ee(:, ii);
    inan = isnan(y0);
    ny = numel(find(~inan));
    if ny > 0
        dat.y(ii) = sum(y0(~inan)) ./ ny;
        dat.e(ii) = sqrt(sum(e0(~inan).^2)) ./ ny;
    end
end