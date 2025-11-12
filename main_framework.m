img_dist = imreadbw('img1.jpg');
imagesc(img_dist)
colormap gray

K = [ 388.6795 0  343.7415;
       0  389.4250 234.6182;
       0 0 1];

KNew = [250 0.0 512;
        0 250 384;
        0.0 0.0 1.0];

function df = radial_distortion_factor_1(r)
    df = (1 / (0.926 * r)) * atan(2*r*tan(0.926/2));
end

function img_undist = undistort(img_dist, K)
    size = size(img_dist)
    H = size(1);
    W = size(2);

    img_undist = zeros(H, W);
    for v = 1:H
        for u = 1:W
            % calculate r
            r = sqrt((v-H)^2 + (u-W)^2);
            df = radial_distortion_factor_1(r);
            p = K * [u; v; 1];
            p_undist = p * df;
            u_undist = round(p_undist(1) / p_undist(3));
            v_undist = round(p_undist(2) / p_undist(3));
            if u_undist > 0 && u_undist <= W && v_undist > 0 && v_undist <= H
                img_undist(v_undist, u_undist) = img_dist(v, u);
            end
        end
    end
end

function img_undist = undistort(img_dist, K, KNew, H_new, W_new, radial_distortion_factor)
    [H_orig, W_orig] = size(img_dist);
    [U_new, V_new] = meshgrid(1:W_new, 1:H_new);

    ones_mat = ones(H_new, W_new);
    P_new = cat(3, U_new, V_new, ones_mat);  % H_new x W_new x 3

    KNew_inv = inv(KNew);
    P_new_reshape = reshape(P_new, [], 3)';  % 3 x (H_new*W_new)
    P_normalized = KNew_inv * P_new_reshape;  % 3 x (H_new*W_new)

    X = P_normalized(1, :);
    Y = P_normalized(2, :);
    r = sqrt(X.^2 + Y.^2);

    df = zeros(size(r));
    for i = 1:length(r)
        df(i) = radial_distortion_factor(r(i));
    end

    X_dist = X .* df;
    Y_dist = Y .* df;
    P_normalized_dist = [X_dist; Y_dist; ones(1, length(X))];

    P_orig = K * P_normalized_dist;  % 3 x (H_new*W_new)

    U_orig = P_orig(1, :) ./ P_orig(3, :);
    V_orig = P_orig(2, :) ./ P_orig(3, :);

    U_orig = reshape(U_orig, H_new, W_new);
    V_orig = reshape(V_orig, H_new, W_new);

    [X_grid, Y_grid] = meshgrid(1:W_orig, 1:H_orig);
    img_undist = interp2(X_grid, Y_grid, double(img_dist), U_orig, V_orig, 'linear', 0);

    img_undist(isnan(img_undist)) = 0;

    if isa(img_dist, 'uint8')
        img_undist = uint8(img_undist);
    end
end

img_undist = undistort(img_dist, K, KNew, 768, 1024, @radial_distortion_factor_1);

subplot(1,2,2)
imagesc(img_dist)
colormap gray
axis equal
subplot(1,2,1)
imagesc(img_undist)
colormap gray
axis equal

%%

img_dist = imreadbw('img2.jpg');

K = [279.7399 0 347.32012;
     0 279.7399 234.99819;
     0 0 1];

KNew = [200 0.0 512;
        0 200 384;
        0.0 0.0 1.0];

function df = radial_distortion_factor_2(r)
    df = 1 - 0.3407 * r + 0.057 * r^2 - 0.0046 * r^3 + 0.00014 * r^4;
end

img_undist = undistort(img_dist, K, KNew, 768, 1024, @radial_distortion_factor_2);
interpn()

subplot(1,2,2)
imagesc(img_dist)
colormap gray
axis equal
subplot(1,2,1)
imagesc(img_undist)
colormap gray
axis equal
