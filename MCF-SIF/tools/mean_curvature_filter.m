function u = mean_curvature_filter(u0, iters, dt)
% Mean-curvature motion smoothing (iterative PDE)
u = u0;
epsv = 1e-6;
for t = 1:iters
    ux = (circshift(u,[0 -1]) - circshift(u,[0 1]))/2;
    uy = (circshift(u,[-1 0]) - circshift(u,[1 0]))/2;
    uxx = circshift(u,[0 -1]) - 2*u + circshift(u,[0 1]);
    uyy = circshift(u,[-1 0]) - 2*u + circshift(u,[1 0]);
    uxy = (circshift(u,[-1 -1]) - circshift(u,[-1 1]) - ...
           circshift(u,[1 -1]) + circshift(u,[1 1]))/4;
    grad2 = ux.^2 + uy.^2;
    denom = (grad2 + epsv).^(3/2);
    kappa = (uxx.*(uy.^2) - 2.*ux.*uy.*uxy + uyy.*(ux.^2)) ./ denom;
    u = u + dt .* kappa .* sqrt(grad2 + epsv);
end
end