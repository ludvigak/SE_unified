function pre = laplace_precomp(opt)

pre = fse_precomp(opt, @kernels.harmonic);