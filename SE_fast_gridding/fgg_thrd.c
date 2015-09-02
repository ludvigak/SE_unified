
static void inline grid_thrd_setup(grid_thrd_ws_t* ws, const int* npdims, int p)
{
#ifdef FGG_THRD
    // Determine which block current thread is working on
    int num_threads = omp_get_num_threads();
    int thread_num = omp_get_thread_num();
    int block_size = npdims[0] / num_threads;
    ws->block_begin = thread_num * block_size;
    ws->block_end = ws->block_begin+block_size;
    if (thread_num+1 == num_threads)
	ws->block_end = npdims[0];
    ws->block_size = npdims[1]*npdims[2];    
    ws->p = p;
#else
    ws->p = p;
#endif
}

static void inline grid_thrd_slice(grid_thrd_ws_t* ws, 
				   int* support_begin, 
				   int* first_slice, 
				   int* last_slice, 
				   int* zidx)
{
#ifdef FGG_THRD
    // Figure out which slices of current P^3 block that go into this block.
    *first_slice = 0;	
    *zidx = 0;
    int block_id = *support_begin / ws->block_size;
    // Skip point if above our chunk
    if (block_id >= ws->block_end)
    {
	ws->skip = true;
	return;
    }
    int diff = ws->block_begin - block_id;
    if (diff > 0)
    {
	// Skip point if too far below our chunk	    
	if (diff > ws->p)
	{
	    ws->skip = true;
	    return;
	}
	// Jump forward to first block inside our chunk
	*support_begin += diff*ws->block_size;
	*zidx += diff*ws->p*ws->p;
	block_id += diff;
	*first_slice += diff;
    }
    // Only iterate over blocks left inside our chunk, but max p of them
    *last_slice = *first_slice + (ws->block_end - block_id);
    if (*last_slice > ws->p)
	*last_slice = ws->p;
    ws->skip = false;
#else
    // Fallback to give normal behavior
    *first_slice = 0;
    *last_slice = ws->p;
    *zidx = 0;
    ws->skip = false;
#endif
}
