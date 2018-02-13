import sys
import matplotlib
import matplotlib.pyplot as plt
import math

def plot(filename, vline):
    with open(filename, 'r') as f:
        line = f.readlines()[3]
        nums = str.split(line, ', ')
        nums = [int(x) for x in nums[:-1]]
        imm = 0.062641
        plt.plot(nums)
        plt.axvline(x=vline, color='r')
        plt.show()

def clip(n, lower, upper):
    return max(lower, min(n, upper))

def attr_bin(attr_value, bin_width, bins):
    return clip(math.floor(attr_value / bin_width), 0, bins)

# args:
# 1: file
# 2: steps
# 3: minval
# 4: maxval
# 5: red_line_loc
steps = pow(2, float(sys.argv[2]))
minval = float(sys.argv[3])
maxval = float(sys.argv[4])
rng = maxval - minval
bin_width = (2 * (0.75*rng - 0.25*rng) / pow(steps, 0.33333)) # TODO 0.75*rng - 0.25*rng looks sus
bins = rng // bin_width


bin_idx = attr_bin(float(sys.argv[5]), bin_width, bins)
plot(sys.argv[1], bin_idx)

exit(0)


