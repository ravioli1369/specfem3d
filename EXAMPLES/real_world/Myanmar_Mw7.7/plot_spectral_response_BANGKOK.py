#!/usr/bin/env python
#
# script to plot figure for spectral responses for Bangkok building collapse site
#
#################################################################################
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # Optional: For more advanced styling

# show figure (or only save to output file)
show = False
add_subplot = True
datadir = "OUTPUT_FILES"  # default input file directory

# get input
for arg in sys.argv:
    if "--dir" in arg:
        datadir = arg.split('=')[1]
    elif "--show" in arg:
        show = True

print("plotting spectral responses:")


# checks data directory
if not os.path.isdir(datadir):
    print(f"  directory {datadir} not found")
    print("")
    print("Please specify directory holding files DB.BANGKOK.BX*.sema.spectral_response.dat with --dir=<my-data-dir>")
    print("exiting...")
    print("")
    sys.exit(1)

files = ( datadir + "/DB.BANGKOK.BXZ.sema.spectral_response.dat",
          datadir + "/DB.BANGKOK.BXN.sema.spectral_response.dat",
          datadir + "/DB.BANGKOK.BXE.sema.spectral_response.dat" )

print("files:")
for file in files:
    print("  ",file)
print("")

dat1 = np.loadtxt(files[0])
dat2 = np.loadtxt(files[1])
dat3 = np.loadtxt(files[2])

x = dat1[:,0]
y1 = dat1[:,1]
y2 = dat2[:,1]
y3 = dat3[:,1]

# --- Enhanced Styling ---

# color palette (using seaborn's default)
sns.set_palette("viridis")  # Try other palettes like "viridis", "muted", "pastel", "deep"
colors = sns.color_palette()  # Get the color palette

# Adjust figure size for better readability
if add_subplot:
    fig, ax1 = plt.subplots(figsize=(10, 7))  # Create the main figure and axes
else:
    plt.figure(figsize=(10, 7))

# colored areas with transparency
alpha_value = 0.3  # Adjust for desired transparency

plt.fill_between(x, y1, 0, color=colors[0], alpha=alpha_value, label='_nolegend_') # _nolegend_ prevents duplicate legend entries
plt.fill_between(x, y2, 0, color=colors[1], alpha=alpha_value, label='_nolegend_')
plt.fill_between(x, y3, 0, color=colors[2], alpha=alpha_value, label='_nolegend_')

# lines
plt.plot(x, y1, label='Z - vertical component', color=sns.color_palette()[0], linestyle='-', linewidth=2, marker='o', markersize=6)
plt.plot(x, y2, label='N - horizontal', color=sns.color_palette()[1], linestyle='-', linewidth=2, marker='s', markersize=6)
plt.plot(x, y3, label='E - horizontal', color=sns.color_palette()[2], linestyle='-', linewidth=2, marker='^', markersize=6)
plt.ylim(0,0.01)

plt.xlabel('Period (s)', fontsize=14)
plt.ylabel('Spectral Acceleration (g)', fontsize=14)
plt.title('Spectral response acceleration - BANGKOK Chatuchak building site', fontsize=16, fontweight='bold')

# Enhance the legend
plt.legend(fontsize=12, loc='upper right', frameon=True, shadow=False)

# Improve grid appearance
plt.grid(True, linestyle='--', alpha=0.5)

# Add annotations (optional, but can highlight specific points)
#plt.annotate('Peak SA', xy=(3.2, 0.006), xytext=(4.5, 0.007),
#             arrowprops=dict(facecolor='black', shrink=0.05), fontsize=10)

# Adjust tick parameters for better readability
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Add a subtle background style (optional)
plt.style.use('ggplot') # Try other styles like 'ggplot', 'seaborn-v0_8-darkgrid', 'seaborn-v0_8-whitegrid'

# --- Create the subplot ---
# Sample data for the extra line plot in the subplot
file = datadir + "/DB.BANGKOK.BXZ.semd"
if add_subplot and os.path.isfile(file):
    # reads trace
    dat = np.loadtxt(file)
    x_subplot = dat[:,0]
    y_subplot = dat[:,1]

    # Define the position and size of the subplot (left, bottom, width, height) as fractions of the figure size
    left, bottom, width, height = 0.6, 0.6, 0.38, 0.2
    ax2 = fig.add_axes([left, bottom, width, height])

    # Show a line at y=0
    #ax2.axhline(y=0, color='black', linewidth=0.5, linestyle='-')

    # Move x-axis ticks to y=0
    ax2.spines['bottom'].set_position(('data', 0))
    ax2.spines['top'].set_visible(False)  # Hide the default top spine
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')

    # Plot a y=0 line in the subplot
    ax2.plot(x_subplot, y_subplot, color=colors[0], linestyle='-', linewidth=1.2, label='displacement - Z component')

    # Style the subplot (optional)
    #ax2.set_xlabel('Time', fontsize=10)
    #ax2.set_ylabel('Value', fontsize=10)
    #ax2.set_title('Additional Data', fontsize=12)
    #ax2.tick_params(axis='both', which='major', labelsize=8)
    ax2.set_xticks([0,500])
    ax2.set_yticks([])
    #ax2.grid(True, linestyle=':', alpha=0.0)
    ax2.grid(False)
    ax2.set_facecolor('white')
    ax2.legend(fontsize=8,loc='lower right')

    # --- Show the plot ---
    plt.tight_layout(rect=[0, 0, 1, 1]) # Adjust layout to make space for the subplot

else:
    # --- Show the plot ---
    plt.tight_layout()  # Adjust layout to prevent labels from overlapping

# saves as JPEG file
filename = datadir + "/spectral_acceleration_BANGKOK.jpg" # "out.jpg"
plt.savefig(filename)
print("")
print("plotted spectral responses in: ",filename)
print("")

if show: plt.show()

