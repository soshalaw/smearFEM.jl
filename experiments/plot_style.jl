# Shared plot styling for the experiment scripts. `include` this from a file that already has
# `using Plots` in scope — the definitions land in the includer's namespace, and `RGB` comes
# from Plots.
#
# `PLOT_CONFIG` deliberately stays per-file: the four scripts that define one do NOT agree
# (font 11 vs 12, 320x480 vs 360x330, and three different margin sets), so a shared dict would
# silently resize their figures. Only the colours are genuinely common.

global def_orange = RGB(245/255,118/255,0)
global def_blue = RGB(5/255,79/255,185/255)
global def_red = RGB(196/255,70/255,1/255)
global def_green = RGB(2/255,147/255,86/255)
