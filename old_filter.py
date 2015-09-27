import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import multiprocessing as mp


file_to_filter = "/media/user/DataFB/AutoHeadFix_Data/0731/EL_LRL/Videos/M1312000377_1438367429.638635.raw"
file_to_save = "/media/user/DataFB/AutoHeadFix_Data/0731/EL_LRL/test.rawf"
corrfile_to_save = "/media/user/DataFB/AutoHeadFix_Data/0731/EL_LRL/corr.rawf"

width = 256
height = 256
frame_rate = 30.0
frame_size = width * height * 3

starting_frame = 100

def get_frames(rgb_file):
    global total_number_of_frames
    with open(rgb_file, "rb") as file:
        frames = np.fromfile(file, dtype=np.uint8)
        total_number_of_frames = int(np.size(frames)/frame_size)

        frames = np.reshape(frames, (total_number_of_frames, width, height, 3))

        frames = frames[starting_frame:, :, :, 1]
        frames = np.asarray(frames, dtype=np.float32)
        total_number_of_frames = frames.shape[0]


    return frames



nyq = frame_rate/2.0
low_limit = 0.1/nyq
high_limit = 1.0/nyq
order = 4
rp = 0.1
Wn = [low_limit, high_limit]

def cheby_filter(frames):
    b, a = signal.cheby1(order, rp, Wn, 'bandpass', analog=False)
    print "Filtering..."
    frames = signal.filtfilt(b, a, frames, axis=0)
    #for i in range(frames.shape[-1]):
    #    frames[:, i] = signal.filtfilt(b, a, frames[:, i])
    print "Done!"
    return frames


def calculate_avg(frames):
    return np.mean(frames, axis=0)


def calculate_df_f0(frames):
    print frames.shape
    baseline = np.mean(frames, axis=0)
    frames = np.divide(np.subtract(frames, baseline), baseline)

    return frames

def save_to_file(filename, frames, dtype):
    with open(filename, "wb") as file:
        frames.astype(dtype).tofile(file)

# Left Barrel
#seed_x = 139
#seed_y = 57

# Left Visual 
#seed_x = 200
#seed_y = 80

# Left Hind Limb
#seed_x = 126
#seed_y = 97

# MC
seed_x = 62
seed_y = 113

##
##def correlation_map(seed_x, seed_y, frames):
##    seed_pixel = np.asarray(frames[:, seed_x, seed_y], dtype=np.float32)
##    
##    print np.shape(seed_pixel)
##    # Reshape into time and space
##    frames = np.reshape(frames, (total_number_of_frames, width*height))
##    print np.shape(frames)
##    correlation_map = []
##    for i in range(frames.shape[-1]):
##        correlation_map.append(pearsonr(frames[:, i], seed_pixel)[0])
##        
##    correlation_map = np.asarray(correlation_map, dtype=np.float32)
##    correlation_map = np.reshape(correlation_map, (width, height))
##    print np.shape(correlation_map)
##
##    return correlation_map
                        

class CorrelationMapDisplayer:
    def __init__(self, frames):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        image = self.get_correlation_map(128, 128, frames)
        
        self.imgplot = self.ax.imshow(image)
        
        self.canvas = self.fig.canvas
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)

        self.frames = frames

        


    def display(self, c_map, low, high):
        self.imgplot.set_cmap(c_map)
        self.imgplot.set_clim(low, high)
        self.fig.colorbar(self.imgplot)

        plt.show()

    def onclick(self, event):
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata)
        # X and Y must be flipped to have correct dimmensions!
        image = self.get_correlation_map(int(event.ydata), int(event.xdata), self.frames)

        self.imgplot = self.ax.imshow(image)
        self.canvas.draw()

    def get_correlation_map(self, seed_x, seed_y, frames):
        seed_pixel = np.asarray(frames[:, seed_x, seed_y], dtype=np.float32)
        
        print np.shape(seed_pixel)
        # Reshape into time and space
        frames = np.reshape(frames, (total_number_of_frames, width*height))
        print np.shape(frames)
        correlation_map = []

        
        for i in range(frames.shape[-1]):
            correlation_map.append(pearsonr(frames[:, i], seed_pixel)[0])
            
        correlation_map = np.asarray(correlation_map, dtype=np.float32)
        correlation_map = np.reshape(correlation_map, (width, height))
        print np.shape(correlation_map)

        return correlation_map



def display_image(image, c_map, low_end_limit, high_end_limit, frames):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    imgplot = ax.imshow(image)
    
    imgplot.set_cmap(c_map)
    imgplot.set_clim(low_end_limit, high_end_limit)
    fig.colorbar(imgplot)

    displayer = CorrelationMapDisplayer(fig, image, frames)

    plt.show()





frames = get_frames(file_to_filter)
avg_pre_filt = calculate_avg(frames)
#frames = np.reshape(frames, (total_number_of_frames, width*height))
frames = cheby_filter(frames)
#frames = np.reshape(frames, (total_number_of_frames, width, height))
frames += avg_pre_filt

frames = calculate_df_f0(frames)

#corr_map = correlation_map(seed_x, seed_y, frames)

mapper = CorrelationMapDisplayer(frames)
mapper.display('spectral', 0.5, 1.0)

#display_image(corr_map, 'spectral', 0.5, 1.0, frames)

#save_to_file(corrfile_to_save, corr_map, np.float32)
#save_to_file(file_to_save, frames, np.float32)


