function pre = stresslet_precomp(opt)

pre = fse_precomp(opt, @kernels.biharmonic);
