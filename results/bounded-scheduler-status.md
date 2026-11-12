# Bounded subdivision scheduler

Temporary interval arrays now contain at most batch_size boxes. A stack retains every remaining chunk and both children of every unresolved box. The arithmetic kernels and the criterion for discarding or certifying a box are unchanged.

This bounds the working batch, not total runtime or all retained memory. The maximum box count remains an explicit failure limit. The seven-charge near-degenerate instance still needs a measured run and independent certificate verification. No certificate is emitted after a budget error.
