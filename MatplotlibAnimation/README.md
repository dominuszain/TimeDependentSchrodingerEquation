# Wave-packet animation

`main.py` animates the probability density in `wavepacket.npz`, which contains
11 sampled points of the wave packet at 10,000 time steps. Run:

```bash
python main.py
```

This writes `output.gif`. Increase the number of sampled points in the Scilab
code for a smoother animation. Long runs use a lot of memory; a practical
approach is to render several short segments and join them afterwards.
