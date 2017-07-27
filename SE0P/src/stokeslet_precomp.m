function pre = stokeslet_precomp(opt)

pre = fse_precomp(opt, @kernels.biharmonic);
