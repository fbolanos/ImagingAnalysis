from filter import *
from math import pow, sqrt
import os
import numpy as np
import matplotlib.pyplot as plt
import cv2

base_dir = ["/media/user/DataFB/AutoHeadFix_Data/0731/EL_LRL/",
            "/media/user/DataFB/AutoHeadFix_Data/0730/EL_LRL_fluc/",
            "/media/user/DataFB/AutoHeadFix_Data/0806/EL_LRL/",
            "/media/user/DataFB/AutoHeadFix_Data/0812/EP_LRL/",
            "/media/user/DataFB/AutoHeadFix_Data/0618/"]
mice = ["1312000377", "1302000245", "1312000300", "2015050202", "2015050115"]

master_filenames = ["/media/user/DataFB/AutoHeadFix_Data/377_master.raw",
                    "/media/user/DataFB/AutoHeadFix_Data/245_master.raw",
                    "/media/user/DataFB/AutoHeadFix_Data/300_master.raw",
                    "/media/user/DataFB/AutoHeadFix_Data/202_master.raw",
                    "/media/user/DataFB/AutoHeadFix_Data/115_master.raw"]



#mice = ["1312000377", "2015050115", "1302000245", "1312000159", "1312000300", "2015050202", "1312000573"]
frame_oi = 400
limit_time = 1438365391.603202

class Position:
    def __init__(self, dx, dy):
        self.dx = dx
        self.dy = dy

        self.dd = sqrt(pow(dx, 2) + pow(dy, 2))

def get_file_list(base_dir, mouse):
    print base_dir
    lof = []
    lofilenames = []
    for root, dirs, files in os.walk(base_dir):
            for file in files:
                if (file.endswith(".raw") or file.endswith(".mp4")) and mouse in file:
                    #print os.path.join(root, file)
                    #if check_time_is_valid(file):
                    #    lof.append((os.path.join(root, file)))
                    #    lofilenames.append(file)
                    #else:
                    #    print file
                    lof.append((os.path.join(root, file)))
                    lofilenames.append(file)
    return lof, lofilenames


def check_time_is_valid(video_file):
    time = video_file[12:-4]
    if float(time) <= limit_time:
        return False
    else:
        return True

def get_video_frames(lof):
    list_of_trial_frames = []
    for video_file in lof:
        print "Getting frames: " + video_file
        frames = get_frames(video_file)
        frames = frames[:800, :, :]
        print np.shape(frames)
        if frames.shape[0] >= 800:
            list_of_trial_frames.append(frames)
        else:
            print "Did not add this file, not enough frames."
    print np.shape(list_of_trial_frames)
    return list_of_trial_frames

    

def get_all_processed_frames(lof):
    list_of_trial_frames = []
    for video_file in lof:
        print "Getting frames: " + video_file
        frames = get_processed_frames(video_file)
        print np.shape(frames)
        list_of_trial_frames.append(frames)
    print np.shape(list_of_trial_frames)
    return list_of_trial_frames

def find_min_ref(lor):
    curr_min = 100
    print np.shape(lor)
    for positions in lor:
        sum = 0
        for position in positions:
            sum += position.dd
            
        print sum
        if curr_min > sum:
            curr_min = sum
            curr_min_positions = positions
            
    print curr_min
    return curr_min_positions

def get_distance_var_mp4(all_frames):
    filtered_frames = []
    print "Filtering the frame of interest for all trials..."
    for frames in all_frames:
        filtered_frames.append(filter2_test(frames, frame_oi))

    print "Getting all the distances.."
    # Get all the distances using all videos as ref point, thus size of matrix is n^2
    list_of_ref = []
    for frame_ref in filtered_frames:
        list_of_positions = []
        res_trials = parmap.map(image_registration.chi2_shift, filtered_frames, frame_ref) 
        # res_trials is array of trials * [dx, dy, edx, edy]
        for res in res_trials:
            list_of_positions.append(Position(res[0], res[1]))
        #for frame in filtered_frames:
        #    dx, dy, edx, edy = image_registration.chi2_shift(frame_ref, frame)
        #    list_of_positions.append(Position(dx, dy))

        list_of_ref.append(list_of_positions)
    print "Finding the min..."
    list_of_positions = find_min_ref(list_of_ref)

    return list_of_positions

def get_distance_var(lof):
    
    list_of_frames = np.asarray(get_video_frames(lof), dtype=np.uint8)
    
    print "Filtering the frame of interest for all trials..."
    # Filter the frame of interest to make vessels obvious.
    filtered_frames = []
    for frames in list_of_frames:
        filtered_frames.append(filter2_test(frames, frame_oi))

    print "Getting all the distances.."
    # Get all the distances using all videos as ref point, thus size of matrix is n^2
    list_of_ref = []
    for frame_ref in filtered_frames:
        list_of_positions = []
        res_trials = parmap.map(image_registration.chi2_shift, filtered_frames, frame_ref) 
        # res_trials is array of trials * [dx, dy, edx, edy]
        for res in res_trials:
            list_of_positions.append(Position(res[0], res[1]))
        #for frame in filtered_frames:
        #    dx, dy, edx, edy = image_registration.chi2_shift(frame_ref, frame)
        #    list_of_positions.append(Position(dx, dy))

        list_of_ref.append(list_of_positions)
    print "Finding the min..."
    list_of_positions = find_min_ref(list_of_ref)

    return list_of_positions

zz
def get_distance_var_from_master_frame(all_frames, master_frames):
    filtered_frames_oi = []
    list_of_positions = []

    print "Getting frame of interest for all frames."
    for frames in all_frames:
        filtered_frames_oi.append(filter2_test(frames, frame_oi))

    master_frame = master_frames[frame_oi]

    print "Getting distances relative to the master frame"
    res_trials = parmap.map(image_registration.chi2_shift, filtered_frames_oi, master_frame)

    for res in res_trials:
        list_of_positions.append(Position(res[0], res[1]))

    return list_of_positions


    
    
class MouseInfo:
    def __init__(self, tag, p2p_x, p2p_y, avg_x, avg_y, n_trials):
        self.tag = tag
        self.p2p_x = p2p_x
        self.p2p_y = p2p_y
        self.avg_x = avg_x
        self.avg_y = avg_y
        self.n_trials = n_trials

def p2p(arr):
    return max(arr)-min(arr)
        
def do_it_all():
    list_mouse_info = []
    for mouse in mice:
        lof, lofilenames = get_file_list(base_dir, mouse)
        print "Lof: ", lof
        lop = get_distance_var(lof)
        dx_trials = []
        dy_trials = []
        for position in lop:
            dx_trials.append(position.dx)
        for position in lop:
            dy_trials.append(position.dy)
        
        peak_x = p2p(dx_trials)
        peak_y = p2p(dy_trials)
        avg_x = np.mean(dx_trials)
        avg_y = np.mean(dy_trials)

        list_mouse_info.append(MouseInfo(mouse, peak_x, peak_y, avg_x, avg_y, len(lop)))
 
    with open(base_dir+"data.tsv", "w") as file:
        file.write("Tag\tp2p_x\tavg_x\tp2p_y\tavg_y\tn_trials\n")
        for info in list_mouse_info:
            file.write(info.tag + "\t" + str(info.p2p_x) + "\t" + str(info.avg_x) + "\t" + str(info.p2p_y) + "\t" + str(info.avg_y) + "\t" + str(info.n_trials) + "\n")

    print "Done it all!"


def process_frames(frames, freq, mouse, dir):
    print "Fixing paint.."
    mouse = mouse[-3:]
    mask = 0
    with open(dir+mouse+"_paint_mask.raw", "rb") as file:
        mask = np.fromfile(file, dtype=np.float32)
        mask = mask.byteswap()
        
    indices = np.squeeze((mask > 0).nonzero())
    paint_frames = np.zeros((frames.shape[0], len(indices)))
    frames = np.reshape(frames, (frames.shape[0], width*height))
    
    for i in range(frames.shape[0]):
        paint_frames[i, :] = frames[i, indices]
    print np.shape(paint_frames)

    mean_paint = np.mean(paint_frames, axis=1)
    mean_paint /= np.mean(mean_paint)

    print np.shape(mean_paint)
    

    frames = np.divide(frames.T, mean_paint)
    frames = frames.T

    frames = np.reshape(frames, (frames.shape[0], width, height))
        
    print "Calculating mean..."
    avg_pre_filt = calculate_avg(frames)

    print "Temporal filter... ", freq.low_limit, "-", freq.high_limit, "Hz"
    frames = cheby_filter(frames, freq.low_limit, freq.high_limit)
    frames += avg_pre_filt

    print "Calculating DF/F0..."
    frames = calculate_df_f0(frames)

    print "Applying MASKED GSR..."
    #frames = gsr(frames)
    frames = masked_gsr(frames, dir+mouse+"_mask.raw")

    #print "Getting SD map..."
    #sd = standard_deviation(frames)

    return frames


def shift_frames(frames, positions):
    print positions.dx, positions.dy
    print frames.shape
    for i in range(len(frames)):
        frames[i] = image_registration.fft_tools.shift2d(frames[i], positions.dx, positions.dy)

    return frames
    

def align_frames(mouse, dir, freq):
    lofiles, lofilenames = get_file_list(dir+"Videos/", mouse)
    print lofilenames
    lop = get_distance_var(lofiles)

    all_frames = np.asarray(get_video_frames(lofiles), dtype=np.uint8)
    print "Alligning all video frames..."

    all_frames = parmap.starmap(shift_frames, zip(all_frames, lop))
##    for i in range(len(lop)):
##        for frame in all_frames[i]:
##            frame = image_registration.fft_tools.shift2d(frame, lop[i].dx, lop[i].dy)

    print np.shape(all_frames)

    count = 0
    new_all_frames = parmap.map(process_frames, all_frames, freq, mouse, dir)
    '''
    for frames in all_frames:
        print np.shape(frames)
        save_to_file("Green/"+lofilenames[count][:-4]+"_aligned.raw", frames, np.float32)

        print "Calculating mean..."
        avg_pre_filt = calculate_avg(frames)

        print "Temporal filter..."
        frames = cheby_filter(frames)
        frames += avg_pre_filt
        save_to_file("Green/Cheby/"+lofilenames[count][:-4]+"_BPFilter_0.1-1Hz.raw", frames, np.float32)


        print "Calculating DF/F0..."
        frames = calculate_df_f0(frames)
        save_to_file("Green/DFF/"+lofilenames[count][:-4]+"_DFF.raw", frames, np.float32)

        print "Applying MASKED GSR..."
        #frames = gsr(frames)
        frames = masked_gsr(frames, save_dir+"202_mask.raw")
        save_to_file("Green/GSR/"+lofilenames[count][:-4]+"_GSR.raw", frames, np.float32)


        print "Getting SD map..."
        sd = standard_deviation(frames)
        save_to_file("Green/SD_maps/"+lofilenames[count][:-4]+"_SD.raw", frames, np.float32)

        new_all_frames.append(frames)
        count += 1
    '''
    print "Creating array..."
    new_all_frames = np.asarray(new_all_frames, dtype=np.float32)
    all_frames = np.asarray(all_frames, dtype=np.float32)
    
    print "Joining Files..."
    new_all_frames = np.reshape(new_all_frames,
                            (new_all_frames.shape[0]*new_all_frames.shape[1],
                            new_all_frames.shape[2],
                            new_all_frames.shape[3]))
    all_frames = np.reshape(all_frames,
                            (all_frames.shape[0]*all_frames.shape[1],
                            all_frames.shape[2],
                            all_frames.shape[3]))

    print "Shapes: "
    print np.shape(all_frames)
    print np.shape(new_all_frames)

    where_are_NaNs = np.isnan(new_all_frames)
    new_all_frames[where_are_NaNs] = 0

    save_to_file("FULL_conc.raw", new_all_frames, np.float32)
    save_to_file("conc_RAW.raw", all_frames, np.float32)
    sd = standard_deviation(new_all_frames)
    save_to_file("FULL_SD.raw", sd, np.float32)

    print "Displaying correlation map..."
    mapper = CorrelationMapDisplayer(new_all_frames)
    mapper.display('spectral', 0.0, 1.0)
    

def process_frames_evoked(frames, freq, mouse, dir):

    print "Fixing paint.."
    mouse = mouse[-3:]
    mask = 0
    with open(dir+mouse+"_paint_mask.raw", "rb") as file:
        mask = np.fromfile(file, dtype=np.float32)
        mask = mask.byteswap()
    
    indices = np.squeeze((mask > 0).nonzero())
    paint_frames = np.zeros((frames.shape[0], len(indices)))
    frames = np.reshape(frames, (frames.shape[0], width*height))
    
    for i in range(frames.shape[0]):
        paint_frames[i, :] = frames[i, indices]
    print np.shape(paint_frames)

    mean_paint = np.mean(paint_frames, axis=1)
    mean_paint /= np.mean(mean_paint)

    print np.shape(mean_paint)
    print "Calculating mean..."
    avg_pre_filt = calculate_avg(frames)

    print "Temporal filter... ", freq.low_limit, "-", freq.high_limit, "Hz"
    frames = cheby_filter(frames, freq.low_limit, freq.high_limit)
    frames += avg_pre_filt

    print "Calculating DF/F0..."
    frames = calculate_df_f0(frames)

    print "Applying MASKED GSR..."
    frames = gsr(frames)
    frames = masked_gsr(frames, dir+mouse+"_mask.raw")


    return frames


def get_evoked_map(mouse, dir, master_frames_filename, freq):
    
    lofiles, lofilenames = get_file_list(dir+"Videos/", mouse)
    print lofilenames
    #all_frames = np.asarray(get_video_frames(lofiles), dtype=np.float32)
    '''
    master_frames = get_frames(master_frames_filename)
    lop = get_distance_var_from_master_frame(all_frames, master_frames)
    print "Alligning all video frames..."
    all_frames = parmap.starmap(shift_frames, zip(all_frames, lop))
    all_frames = np.asarray(all_frames, dtype=np.float32)
    print np.shape(all_frames)
    '''

    lop = get_distance_var(lofiles)
    all_frames = np.asarray(get_video_frames(lofiles), dtype=np.float32)
    print "Alligning all video frames..."
    all_frames = parmap.starmap(shift_frames, zip(all_frames, lop))
    print np.shape(all_frames)
    all_frames = np.asarray(all_frames, dtype=np.float32)

    #all_frames = np.reshape(all_frames,
    #                        (all_frames.shape[0]*all_frames.shape[1],
    #                        all_frames.shape[2],
    #                        all_frames.shape[3])) 
    #save_to_file(dir,"raw_conc_"+mouse+".raw", all_frames, np.float32)


    new_all_frames = parmap.map(process_frames_evoked, all_frames, freq, mouse, dir)
    
    all_frames = np.reshape(all_frames,
                            (all_frames.shape[0]*all_frames.shape[1],
                            all_frames.shape[2],
                            all_frames.shape[3])) 


    print "Creating array.."
    new_all_frames = np.asarray(new_all_frames, dtype=np.float32)
    print "Averaging together..."
    new_all_frames = np.mean(new_all_frames, axis=0)

    print np.shape(new_all_frames)
    #count = 1
    #for frames in new_all_frames:
    #    print "Saving trial %s." % (str(count))
    #    save_to_file("/media/user/DataFB/AutoHeadFix_Data/0815/EL_noled/stim/", "trial_"+str(count)+".raw", frames, np.float32)
    #    count += 1
    
    #save_to_file(dir,"raw_conc_"+mouse+".raw", all_frames, np.float32)
    save_to_file(dir,"evoked_"+mouse+"_.raw", new_all_frames, np.float32)
    sd = standard_deviation(new_all_frames)
    save_to_file(dir,"all_frames_SD"+mouse+"_.raw", sd, np.float32)


def get_mp4_frames(filename):
    list_of_frames = []
    print filename
    vidcap = cv2.VideoCapture(filename)
    frames = []
    while True:
        success, image = vidcap.read()
        if not success:
            break
        image = np.asarray(image)
        image = image[:, :, 1]
        frames.append(image)

    frames = np.asarray(frames, dtype=np.float32)
    frames = frames[10:710, :, :]
    print np.shape(frames)

    return frames


def process_mp4_frames(frames, freq, mouse, dir):
    mouse = mouse[-3:]
    print "Calculating mean..."
    avg_pre_filt = calculate_avg(frames)

    print "Temporal filter... ", freq.low_limit, "-", freq.high_limit, "Hz"
    frames = cheby_filter(frames, freq.low_limit, freq.high_limit)
    frames += avg_pre_filt

    print "Calculating DF/F0..."
    frames = calculate_df_f0(frames)

    print "Applying MASKED GSR..."
    frames = masked_gsr(frames, dir+mouse+"_mask.raw")


    return frames

def get_corr_maps_mp4(mouse, dir, freq):
    str_freq = str(freq.low_limit) + "-" + str(freq.high_limit) + "Hz"
    lofiles, lofilenames = get_file_list(dir+"MP4/", mouse)
    print lofilenames
    all_frames = []
    for filename in lofiles:
        print "Opening " + filename
        all_frames.append(get_mp4_frames(filename))

    all_frames = np.asarray(all_frames, dtype=np.float32)
    lop = get_distance_var_mp4(all_frames)

    print "Alligning all video frames..."
    all_frames = parmap.starmap(shift_frames, zip(all_frames, lop))
    all_frames = np.asarray(all_frames, dtype=np.float32)
    print np.shape(all_frames)
    print "Joining Files..."
    number_of_files = all_frames.shape[0]
    frames_per_file = all_frames.shape[1]
    all_frames = np.reshape(all_frames,
                            (number_of_files*frames_per_file,
                            width,
                            height))
    print "Saving raw concatanated frames..."
    save_to_file(dir, "115.raw", all_frames, np.float32)
    all_frames = np.reshape(all_frames,
                        (number_of_files,
                         frames_per_file,
                        width,
                        height))

    
    all_frames = parmap.map(process_mp4_frames, all_frames, freq, mouse, dir)
    all_frames = np.asarray(all_frames, dtype=np.float32)

    all_frames = np.reshape(all_frames,
                            (number_of_files*frames_per_file,
                            width,
                            height))
    print "Saving processed frames..."
    save_to_file(dir, "115_processed.raw", all_frames, np.float32)

    

    print "Displaying correlation map..."
    mapper = CorrelationMapDisplayer(all_frames, dir, mouse)
    mapper.display('spectral', 0.0, 1.0)

    print "All done!! :))"


def get_corr_maps(mouse, dir, freq, coords, master_frames_filename):
    str_freq = str(freq.low_limit) + "-" + str(freq.high_limit) + "Hz"
    lofiles, lofilenames = get_file_list(dir+"Videos/", mouse)
    print lofilenames

    #all_frames = np.asarray(get_video_frames(lofiles), dtype=np.float32)
    #master_frames = get_frames(master_frames_filename)
    #lop = get_distance_var_from_master_frame(all_frames, master_frames)
    #print "Alligning all video frames..."
    #all_frames = parmap.starmap(shift_frames, zip(all_frames, lop))


    lop = get_distance_var(lofiles)
    all_frames = np.asarray(get_video_frames(lofiles), dtype=np.uint8)
    print "Alligning all video frames..."
    all_frames = parmap.starmap(shift_frames, zip(all_frames, lop))
    print np.shape(all_frames)

    count = 0
    new_all_frames = parmap.map(process_frames, all_frames, freq, mouse, dir)

    print "Creating array..."
    new_all_frames = np.asarray(new_all_frames, dtype=np.float32)
    all_frames = np.asarray(all_frames, dtype=np.float32)
    
    print "Joining Files..."
    new_all_frames = np.reshape(new_all_frames,
                            (new_all_frames.shape[0]*new_all_frames.shape[1],
                            new_all_frames.shape[2],
                            new_all_frames.shape[3]))

    print "Shapes: "
    print np.shape(all_frames)
    print np.shape(new_all_frames)

    where_are_NaNs = np.isnan(new_all_frames)
    new_all_frames[where_are_NaNs] = 0
    
    save_to_file(dir,"raw_conc_"+mouse+"_"+str_freq+".raw", all_frames, np.float32)
    print "Saving the processed concatenated file..."
    save_to_file(dir,"processed_conc_"+mouse+"_"+str_freq+".raw", new_all_frames, np.float32)
    #sd = standard_deviation(new_all_frames)
    #save_to_file(dir,"all_frames_SD"+mouse+"_"+str_freq+".raw", sd, np.float32)

    #for coord in coords:
    #    corr_map = get_correlation_map(coord.x, coord.y, new_all_frames)
    #    save_to_file(dir, "All_Maps/"+mouse+"_map_"+str(coord.x)+","+str(coord.y)+"_"+str_freq+".raw", corr_map, dtype=np.float32)

    print "Displaying correlation map..."
    mapper = CorrelationMapDisplayer(new_all_frames, dir, mouse)
    mapper.display('spectral', 0.0, 1.0)

    print "All done!! :))"


class FrequencyLimit:
    def __init__(self, low, high):
        self.low_limit = low
        self.high_limit = high


class Coordinate:
    def __init__(self, x, y):
        self.x = x
        self.y = y



def get_correlation_map(seed_x, seed_y, frames):
    seed_pixel = np.asarray(frames[:, seed_x, seed_y], dtype=np.float32)
    
    print np.shape(seed_pixel)
    # Reshape into time and space
    frames = np.reshape(frames, (frames.shape[0], width*height))
    print np.shape(frames)
    print 'Getting correlation... x=', seed_x, ", y=", seed_y

    correlation_map = parmap.map(corr, frames.T, seed_pixel)
    correlation_map = np.asarray(correlation_map, dtype=np.float32)
    correlation_map = np.reshape(correlation_map, (width, height))
    print np.shape(correlation_map)

    return correlation_map

    
#do_it_all()
frequencies = [FrequencyLimit(0.3, 3.),
               FrequencyLimit(0.01, 6.0)]
coords = [Coordinate(138, 192)]

##for i in range(3, len(mice)):
##     for freq in frequencies:
##         get_corr_maps(mice[i], base_dir[i], freq, coords)

get_corr_maps_mp4(mice[4], base_dir[4], frequencies[0])

#get_corr_maps(mice[4], base_dir[4], frequencies[0], coords, master_filenames[3])
#get_evoked_map(mice[4], base_dir[4], master_filenames[3], frequencies[1])


#freq = FrequencyLimit(0.3, 3.0)
#align_frames(mice[0], base_dir[0], freq)


#get_evoked_map(mice[6])

##test_gcamp = get_frames("/media/user/DataFB/AutoHeadFix_Data/0731/EL_LRL/Videos/M1312000377_1438367187.563086.raw")
##
##avg_pre_filt = calculate_avg(test_gcamp)
##test_gcamp = cheby_filter(test_gcamp)
##test_gcamp += avg_pre_filt
##test_gcamp = calculate_df_f0(test_gcamp)
##test_gcamp = masked_gsr(test_gcamp, save_dir+"377_mask.raw")


#test_gfp = get_frames("/media/user/DataFB/AutoHeadFix_Data/0806/EL_LRL/Videos/M1312000300_1438887348.410214.raw")

#avg_pre_filt = calculate_avg(test_gfp)
#test_gfp = cheby_filter(test_gfp)
#test_gfp += avg_pre_filt
#test_gfp = calculate_df_f0(test_gfp)
#test_gfp = masked_gsr(test_gfp, "/media/user/DataFB/AutoHeadFix_Data/0806/EL_LRL/300_mask.raw")


#mapper= CorrelationMapDisplayer(test_gcamp)
#mapper.display('spectral', -0.5, 1.0)

#mapper2= CorrelationMapDisplayer(test_gfp)
#mapper2.display('spectral', -0.5, 1.0)
                        

#lofiles, lofilenames = get_file_list("/media/user/DataFB/AutoHeadFix_Data/0731/EL_LRL/Green/GSR/", mice[0])
#all_frames = np.asarray(get_all_processed_frames(lofiles), dtype=np.float32)
#print "Joining Files..."
#all_frames = np.reshape(all_frames,
#                        (all_frames.shape[0]*all_frames.shape[1],
#                        all_frames.shape[2],
#                        all_frames.shape[3]))
#print np.shape(all_frames)


#print "Displaying correlation map..."
#mapper = CorrelationMapDisplayer(all_frames)
#mapper.display('spectral', -0.5, 1.0)
#lof = get_file_list(base_dir, mice[1])

# List of positions for all trials with the best reference point
#lop =get_distance_var(lof)

#dx_trials = []
#for position in lop:
#    dx_trials.append(position.dx)
#plt.plot(dx_trials)
#plt.title("Change in X for all trials")
#plt.ylabel("dx in pixels")
#plt.xlabel("Trial number")
#plt.show()

#dy_trials = []
#for position in lop:
#    dy_trials.append(position.dy)
#plt.plot(dy_trials)
#plt.title("Change in Y for all trials")
#plt.ylabel("dy in pixels")
#plt.xlabel("Trial number")
#plt.show()
