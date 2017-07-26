function times = time_min(times, t)
times.pre   = min(t.pre, times.pre);
times.fft   = min(t.fft, times.fft);
times.grid  = min(t.grid, times.grid);
times.int   = min(t.int, times.int);
times.scale = min(t.scale, times.scale);