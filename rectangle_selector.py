from __future__ import print_function
"""
Do a mouseclick somewhere, move the mouse to some destination, release
the button.  This class gives click- and release-events and also draws
a line or a box from the click-point to the actual mouseposition
(within the same axes) until the button is released.  Within the
method 'self.ignore()' it is checked wether the button from eventpress
and eventrelease are the same.

"""
from matplotlib.widgets import RectangleSelector
import numpy as np
from mpl_toolkits.axes_grid.axislines import Subplot
import matplotlib.pyplot as plt
import matplotlib
from filter import *
import os
import scipy

font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 25}

matplotlib.rc('font', **font)
#matplotlib.rcParams.update({'font.size': 18})

list_df_f0 = []
all_true_means = []
all_list_sems = []

foi = "/media/user/DataFB/AutoHeadFix_Data/0815/EL_noled/subset_pros_300.raw"
folder  = "/media/user/DataFB/AutoHeadFix_Data/0815/EL_noled/stim_300/"

width = 256
height = 256
no_frames = 800
start_frame = 495
end_frame = 550

timings = [1438277955.2506770,
            1438277960.400527,
            1438277965.499913,
            1438277970.58546,
            1438277975.708472,
            1438277980.777134]

headfixing_time = timings[0]
for i in range(len(timings)):
    timings[i] = timings[i] - headfixing_time


def get_trial_list(base_dir):
    print (base_dir)
    lof = []
    lofilenames = []
    for root, dirs, files in os.walk(base_dir):
            for file in files:
                lof.append((os.path.join(root, file)))
                lofilenames.append(file)
    return lof, lofilenames


def get_trial_frames(filename):
    with open(filename, "rb") as file:
        frames = np.fromfile(file, dtype=np.float32)
        frames = np.reshape(frames, (no_frames, width, height))
        frames = np.asarray(frames, dtype=np.float32)
        frames = frames[start_frame:end_frame, :, :]

        print (np.shape(frames))

    return frames


def line_select_callback(eclick, erelease):
    global x1, x2, y1, y2
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print ("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    print (" The button you used were: %s %s" % (eclick.button, erelease.button))

    

def toggle_selector(event):
    print (' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print (' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print (' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)
    if event.key in ['T', 't']:
        add_trials_to_avg()
        #add_to_average()
        print ("Added to the trace!")

def add_trials_to_avg():
    global all_true_means, all_list_sems

    list_means = [] # has shape (trial, frames)
    for frames in all_frames:
        means_roi = [] # has shape (frames)
        for frame in frames:
            mean_width = np.mean(frame[y1:y2, x1:x2], axis=0)
            mean_height = 100 * np.mean(mean_width)
            means_roi.append(mean_height)

        list_means.append(means_roi)

    list_means = np.asarray(list_means)

    list_sems = []
    for means_frame in list_means.T:
        list_sems.append(scipy.stats.sem(means_frame))


    true_means = np.mean(list_means, axis=0)
    print (np.shape(true_means))
    print (np.shape(list_sems))
    all_true_means.append(true_means)
    all_list_sems.append(list_sems)
    

        

def add_to_average():
    global list_df_f0

    df_f0 = []
    for frame in frames:
        mean_width = np.mean(frame[y1:y2, x1:x2], axis=0)
        mean_height = np.mean(mean_width)
        df_f0.append(mean_height)

    list_df_f0.append(df_f0)


def plot_test():
    for df_df0 in list_df_f0:
        plt.plot(df_df0)

    plt.show()
    
def plot_df_f0s():
    fig = plt.figure(1, (3,3))

    ax = Subplot(fig, 111)
    fig.add_subplot(ax)

    ax.axis["right"].set_visible(False)
    ax.axis["top"].set_visible(False)

    
    x = np.linspace(-333.3, 1498.5, num=55)
    #x = range(1, 56)
    for true_means, list_sems in zip(all_true_means, all_list_sems):
        ax.errorbar(x, true_means, yerr=list_sems, elinewidth=3, linewidth=4)
    ax.axvspan(0, 400, color='y', alpha=0.5, lw=0)
    ax.set_xticks([-200, 0, 200, 400, 600, 800, 1000, 1200, 1400])
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(r"$\Delta F/F_0$(%)")
    ax.xaxis.set_tick_params(width=2, length=7)
    ax.yaxis.set_tick_params(width=2, length=7)
    
    #x_axis = np.linspace(3.333, 26.6667, num=800) 
    #for df_df0 in list_df_f0:
        #plt.plot(df_df0)
    ax.set_ylim([-1.5,2.5])
    ax.set_xlim([-340, 1500])

    # axes color.... jeff didn't like it :(
    #plt.axhline(linewidth=3, color='k')
    #plt.axvline(linewidth=3, color='k')

    
    #plt.xlim([3.333,26.6667])

    #plt.plot([83, 83], [-0.015, 0.025], 'k-', lw=2)
    #plt.plot([203, 203], [-0.015, 0.025], 'k-', lw=2)
    #plt.plot([505, 505], [-0.015, 0.025], 'k-', lw=2)
    #plt.axvspan(5, 17, color='y', alpha=0.5, lw=0)
    #plt.axvspan(505, 517, color='y', alpha=0.5, lw=0)

    #for i in range(1, len(timings)):
        #plt.plot([timings[i], timings[i]], [-0.015, 0.025], 'k-', lw=2)
    plt.show()
        



lof, lofilnes = get_trial_list(folder)

all_frames = []
for fd in lof:
    all_frames.append(get_trial_frames(fd))

print (np.shape(all_frames))


#frames = get_processed_frames(foi)
#frames = frames.byteswap()
#where_are_NaNs = np.isnan(frames)
#frames[where_are_NaNs] = 0

#for frame in frames:
#    print (np.max(frame) - np.min(frame))

fig, current_ax = plt.subplots()                    # make a new plotingrange
plt.imshow(np.mean(all_frames, axis=0)[35], vmin=-0.02, vmax=0.02)


print ("\n      click  -->  release")

# drawtype is 'box' or 'line' or 'none'
toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                       drawtype='box', useblit=True,
                                       button=[1,3], # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels')
plt.connect('key_press_event', toggle_selector)
plt.show()

plot_df_f0s()
#plot_test()
