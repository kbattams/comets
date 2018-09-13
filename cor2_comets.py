import numpy as np
from sunpy.map import Map
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import glob
from scipy.ndimage.interpolation import rotate
import os
import astropy.units as u

#import cv2

head = [None,None]
tail = [None,None]

def line_select_callback(eclick, erelease):
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    head[:] = x1, y1
    tail[:] = x2, y2
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    print(" The button you used were: %s %s" % (eclick.button, erelease.button))


def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)

def get_rotation(point1, point2):
    dx = np.abs(point1[0]-point2[0])
    dy = np.abs(point1[1]-point2[1])
    theta = np.arctan(dy/dx)
    return np.degrees(theta)+90.

# Adaptive Box Photometry Class
class ABPhotomPt:
    def __init__(self, ObjDatetime, target, backgrounds, target_meta, background_meta):
        self.cometIm = target
        self.bkgIm = backgrounds
        self.cometMeta = target_meta
        self.bkgMeta = background_meta
        self.ID = ObjDatetime           # datatime object for target observation
    description="Object to store comet photometry data"
    author="Karl Battams (NRL)"
    
    def date(self):
        return self.ObjDatetime.strftime('%Y-%m-%d %H:%M:%S')
        
    def show(self):
        plt.imshow(self.target)
        plt.show(block=True)
        
    def getBkg(self):
        return np.median(self.backgrounds, axis=0)

# SETUP
datadir = '/Users/battams/Work/COMETS/PROJECTS/COR2_Comets/PYTHON/testdata/'
box_width = 9           # Default comet width in pixels
half_box = np.floor(box_width/2).astype('int')
num_bkgs = 6            # how many backgrounds to include (excluding self)
bkg_buf = num_bkgs/2    # how many bkgs on each side of current

# Read all files (assume calibrated)
# Get file list
flist = glob.glob(os.path.join(datadir,'*.fts'))

nf = len(flist)


print('Found %i files.' % nf)

# Create an empty data cubes as sunpy Map 
# (ref http://docs.sunpy.org/en/stable/code_ref/map.html)
alldata = Map(flist, sequence=True)
allhdrs = alldata.all_meta()

imsize = alldata
# MAY NEED PRE-PROCESSING HERE FOR EASE OF IDENTIFYING COMET

im=alldata[3].data
#r = cv2.selectROI(im)

fig, current_ax = plt.subplots()                 # make a new plotting range

plt.imshow(im,cmap=plt.get_cmap('jet'),origin='lower')

print("Select bounding box from head to tail")

print("\n      click  -->  release")

# drawtype is 'box' or 'line' or 'none'
toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                       drawtype='line', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True)
plt.connect('key_press_event', toggle_selector)
plt.show(block=True)  # The block command prevents code from continuing until window closed


# Get center of sun in pixel coords (NOT USED - JUST FYI FOR ACCESSING METADATA)
suncen=[None, None]
suncen[:] = alldata[0].reference_pixel[0].value, alldata[0].reference_pixel[1].value

#Calculate necessary rotation angle to make comet vertical
rotang = get_rotation(head, tail)
head_int = np.asarray(head).astype('int')
tail_int = np.asarray(tail).astype('int')

no_crop=0 # cropping flag
if ( np.abs(tail_int[0] - head_int[0]) < box_width ):
    tail_int[0]+=half_box
    head_int[0]-=half_box
    no_crop=1                       # Crop shouldn't be needed
 
if ( np.abs(tail_int[1] - head_int[1]) < box_width ):
    tail_int[1]+=half_box
    head_int[1]-=half_box
    no_crop=1                       # Crop shouldn't be needed
    
# Extract the sub image of just the comet
sub_im = im[head_int[1]:tail_int[1],head_int[0]:tail_int[0]]

# Rotate the sub_img, making sure to pad the image so rotation doesn't cut corners
rot_im = rotate(sub_im, rotang, reshape=True)

# Crop comet down to some width, ensuring we keep indices as INTEGERS
center_pix = np.int(rot_im.shape[1] / 2)

if ~no_crop:
    crop_box = rot_im[:,center_pix-half_box+1:center_pix+half_box+1]

# Apply same transform to all images in periphery
# **** DO THIS *****
curr = 3  # current index; make this generic
ind_list = np.arange(nf)
bkg_inds = np.array([curr-3,curr-2,curr-1,curr+1,curr+2,curr+3])

backgrounds = np.empty((num_bkgs, rot_im.shape[0], rot_im.shape[1]))
subs_meta = [alldata[i].meta for i in bkg_inds]
bkg_images = [ rotate(alldata[i].data[head_int[1]:tail_int[1],head_int[0]:tail_int[0]], rotang, reshape=True) for i in bkg_inds]

#for i in range(num_bkgs):
#    tmp = alldata[bkg_inds[i]].data
#    meta_subs[i] = alldata[bkg_inds[i]].meta
#    tmp_sub = tmp[head_int[1]:tail_int[1],head_int[0]:tail_int[0]]
#    backgrounds[i,:,:] = rotate(tmp_sub, rotang, reshape=True)

# Prepare data for adding to class
com_meta = alldata[3].meta
date_time = alldata[3].date

target_image = rot_im

photPt = ABPhotomPt(date_time, target_image, bkg_images, com_meta, subs_meta)





    
    
    
    
