import sys
import matplotlib
import matplotlib.pyplot as plt

# Instructions:
# first argument is file that contains histogram output, as it is printed
# out by chain in commit c570e24628ab3dcc039574aadf90002206c105f0. Second
# argument is x-position of the red line to be drawn.

def plot(filename, vline):
    with open(filename, 'r') as f:
        line = f.readlines()[3]
        nums = str.split(line, ', ')
        nums = [int(x) for x in nums[:-1]]
        imm = 0.062641
        plt.plot(nums)
        plt.axvline(x=vline, color='r')
        plt.show()

plot(sys.argv[1], float(sys.argv[2]))

exit(0)


