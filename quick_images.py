#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
import numpy as np
from astropy.visualization import simple_norm
import pylab as plt
from astropy import wcs
import os
from reproject import reproject_interp
import reproject
import PIL
#import pyavm
import shutil


# In[2]:


def save_rgb(img, filename, avm=None, flip=-1):
    img = (img*256)
    img[img<0] = 0
    img[img>255] = 255
    img = img.astype('uint8')
    img = PIL.Image.fromarray(img[::flip,:,:])
    img.save(filename)
    if avm is not None:
        base = os.path.basename(filename)
        dir = os.path.dirname(filename)
        avmname = os.path.join(dir, 'avm_'+base)
        avm.embed(filename, avmname)
        shutil.move(avmname, filename)
    return img


# In[17]:


image_filenames_pipe ={
    "f150w": "/orange/adamginsburg/jwst/sgrb2/NB/F150W/pipeline/jw05365-o001_t001_nircam_clear-f150w-merged_i2d.fits",
    "f182m": "/orange/adamginsburg/jwst/sgrb2/NB/F182M/pipeline/jw05365-o001_t001_nircam_clear-f182m-merged_i2d.fits",
    "f187n": "/orange/adamginsburg/jwst/sgrb2/NB/F187N/pipeline/jw05365-o001_t001_nircam_clear-f187n-merged_i2d.fits",
    "f210m": "/orange/adamginsburg/jwst/sgrb2/NB/F210M/pipeline/jw05365-o001_t001_nircam_clear-f210m-merged_i2d.fits",
    "f212n": "/orange/adamginsburg/jwst/sgrb2/NB/F212N/pipeline/jw05365-o001_t001_nircam_clear-f212n-merged_i2d.fits",
    "f300m": "/orange/adamginsburg/jwst/sgrb2/NB/F300M/pipeline/jw05365-o001_t001_nircam_clear-f300m-merged_i2d.fits",
    "f360m": "/orange/adamginsburg/jwst/sgrb2/NB/F360M/pipeline/jw05365-o001_t001_nircam_clear-f360m-merged_i2d.fits",
    "f405n": "/orange/adamginsburg/jwst/sgrb2/NB/F405N/pipeline/jw05365-o001_t001_nircam_clear-f405n-merged_i2d.fits",
    "f410m": "/orange/adamginsburg/jwst/sgrb2/NB/F410M/pipeline/jw05365-o001_t001_nircam_clear-f410m-merged_i2d.fits",
    "f466n": "/orange/adamginsburg/jwst/sgrb2/NB/F466N/pipeline/jw05365-o001_t001_nircam_clear-f466n-merged_i2d.fits", # weird, the filename is different from what is downloaded with the STScI pipeline...
    "f480m": "/orange/adamginsburg/jwst/sgrb2/NB/F480M/pipeline/jw05365-o001_t001_nircam_clear-f480m-merged_i2d.fits",
    "f770w": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o002_t002_miri_f770w/jw05365-o002_t002_miri_f770w_i2d.fits",
    "f1280w": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o002_t002_miri_f1280w/jw05365-o002_t002_miri_f1280w_i2d.fits",
    "f2550w": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o002_t002_miri_f2550w/jw05365-o002_t002_miri_f2550w_i2d.fits",
    #"f770w": "/orange/adamginsburg/jwst/sgrb2/NB/data_reprojected/jw05365-o002_t002_miri_f770w_i2d_pipeline_v0.1_reprj_f466.fits",
    #"f1280w": "/orange/adamginsburg/jwst/sgrb2/NB/data_reprojected/jw05365-o002_t002_miri_f1280w_i2d_pipeline_v0.1_reprj_f466.fits",
    #"f2550w": "/orange/adamginsburg/jwst/sgrb2/NB/data_reprojected/jw05365-o002_t002_miri_f2550w_i2d_pipeline_v0.1_reprj_f466.fits",
}

image_sub_filenames_pipe = {
    "f405n-f410m": "/orange/adamginsburg/jwst/sgrb2/NB/F405_minus_F410cont_pipeline_v0.1.fits",
    "f410m-f405n": "/orange/adamginsburg/jwst/sgrb2/NB/F410_minus_F405_fractional_bandwidth_pipeline_v0.1.fits",
    "f212n-f210m": "/orange/adamginsburg/jwst/sgrb2/NB/F212_minus_F210cont_pipeline_v0.1.fits",
    "f187n-f182m": "/orange/adamginsburg/jwst/sgrb2/NB/F187_minus_F182cont.fits",

}

new_basepath = '/orange/adamginsburg/jwst/sgrb2/NB/data_reprojected/'
repr466_image_filenames = {x: y.replace("i2d", "i2d_pipeline_v0.1_reprj_f466") for x,y in image_filenames_pipe.items()}
repr466_image_filenames = {x: (new_basepath+os.path.basename(y)) for x,y in repr466_image_filenames.items()}
repr466_image_sub_filenames = {x: y.replace("i2d", "i2d_reprj_f466") for x,y in image_sub_filenames_pipe.items()}
repr466_image_sub_filenames = {x: (new_basepath+os.path.basename(y)) for x,y in repr466_image_sub_filenames.items()}


# In[19]:


tgt_header = fits.getheader(image_filenames_pipe['f466n'], ext=('SCI', 1))
#tgt_header = fits.getheader('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f480m/jw05365-o001_t001_nircam_clear-f480m_i2d.fits', ext=('SCI', 1))
#wcs.WCS(tgt_header)
#AVM = pyavm.AVM.from_header(tgt_header)


# In[20]:


for filtername in image_sub_filenames_pipe:
    if not os.path.exists(repr466_image_sub_filenames[filtername]):
        print(f"Reprojecting {filtername} {image_sub_filenames_pipe[filtername]} to {repr466_image_sub_filenames[filtername]}")
        result,_ = reproject_interp(image_sub_filenames_pipe[filtername], tgt_header, hdu_in='SCI')
        hdu = fits.PrimaryHDU(data=result, header=tgt_header)
        hdu.writeto(repr466_image_sub_filenames[filtername], overwrite=True)


# In[21]:


for filtername in image_filenames_pipe:
    if not os.path.exists(repr466_image_filenames[filtername]):
        print(f"Reprojecting {filtername} {image_filenames_pipe[filtername]} to {repr466_image_filenames[filtername]}")
        result,_ = reproject_interp(image_filenames_pipe[filtername], tgt_header, hdu_in='SCI')
        hdu = fits.PrimaryHDU(data=result, header=tgt_header)
        hdu.writeto(repr466_image_filenames[filtername], overwrite=True)


# In[ ]:





# In[7]:


image_filenames ={
    "f150w": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f150w/jw05365-o001_t001_nircam_clear-f150w_i2d.fits",
    "f182m": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f182m/jw05365-o001_t001_nircam_clear-f182m_i2d.fits",
    "f187n": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f187n/jw05365-o001_t001_nircam_clear-f187n_i2d.fits",
    "f210m": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f210m/jw05365-o001_t001_nircam_clear-f210m_i2d.fits",
    "f212n": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f212n/jw05365-o001_t001_nircam_clear-f212n_i2d.fits",
    "f300m": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f300m/jw05365-o001_t001_nircam_clear-f300m_i2d.fits",
    "f360m": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f360m/jw05365-o001_t001_nircam_clear-f360m_i2d.fits",
    "f405n": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_f405n-f444w/jw05365-o001_t001_nircam_f405n-f444w_i2d.fits",
    "f410m": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f410m/jw05365-o001_t001_nircam_clear-f410m_i2d.fits",
    "f466n": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_f444w-f466n/jw05365-o001_t001_nircam_f444w-f466n_i2d.fits",
    "f480m": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f480m/jw05365-o001_t001_nircam_clear-f480m_i2d.fits",
    "f770w": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o002_t002_miri_f770w/jw05365-o002_t002_miri_f770w_i2d.fits",
    "f1280w": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o002_t002_miri_f1280w/jw05365-o002_t002_miri_f1280w_i2d.fits",
    "f2550w": "/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o002_t002_miri_f2550w/jw05365-o002_t002_miri_f2550w_i2d.fits",
}
image_sub_filenames = {
    "f182m-f187n": "/orange/adamginsburg/jwst/sgrb2/NB/F182_minus_F187_fractional_bandwidth.fits",
    "f187n-f182m": "/orange/adamginsburg/jwst/sgrb2/NB/F187_minus_F182cont.fits",
    "f210m-f212n": "/orange/adamginsburg/jwst/sgrb2/NB/F210_minus_F212_fractional_bandwidth.fits",
    "f212n-f210m": "/orange/adamginsburg/jwst/sgrb2/NB/F212_minus_F210cont.fits",
    "f405n-f410m": "/orange/adamginsburg/jwst/sgrb2/NB/F405_minus_F410cont.fits",
    "f410m-f405n": "/orange/adamginsburg/jwst/sgrb2/NB/F410_minus_F405_fractional_bandwidth.fits",
    # TODO "f466n-f410m": "/orange/adamginsburg/jwst/sgrb2/NB/F466_minus_F410.fits",
    "f466n-f405n": "/orange/adamginsburg/jwst/sgrb2/NB/F466_minus_F405.fits",
    "f466n-f480m": "/orange/adamginsburg/jwst/sgrb2/NB/F466_minus_F480cont.fits",
    "f480m-f466n": "/orange/adamginsburg/jwst/sgrb2/NB/F480_minus_F466_fractional_bandwidth.fits",
    "f480m-f300m": "/orange/adamginsburg/jwst/sgrb2/NB/F480-F300.fits",
    "f480m-f360m": "/orange/adamginsburg/jwst/sgrb2/NB/F480M-F360M_no_bandpass_correction.fits",
}
new_basepath = '/orange/adamginsburg/jwst/sgrb2/NB/data_reprojected/'
repr480_image_filenames = {x: y.replace("i2d", "i2d_reprj_f480") for x,y in image_filenames.items()}
repr480_image_filenames = {x: (new_basepath+os.path.basename(y)) for x,y in repr480_image_filenames.items()}
repr480_image_sub_filenames = {x: y.replace("i2d", "i2d_reprj_f480") for x,y in image_sub_filenames.items()}
repr480_image_sub_filenames = {x: (new_basepath+os.path.basename(y)) for x,y in repr480_image_sub_filenames.items()}


# In[10]:


import pyavm


# In[ ]:


tgt_header = fits.getheader(image_filenames['f480m'], ext=('SCI', 1))
tgt_header = fits.getheader('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f480m/jw05365-o001_t001_nircam_clear-f480m_i2d.fits', ext=('SCI', 1))
#wcs.WCS(tgt_header)
AVM = pyavm.AVM.from_header(tgt_header)


# In[12]:


for filtername in image_filenames:
    if not os.path.exists(repr480_image_filenames[filtername]):
        print(f"Reprojecting {filtername} {image_filenames[filtername]} to {repr480_image_filenames[filtername]}")
        result,_ = reproject.reproject_interp(image_filenames[filtername], tgt_header, hdu_in='SCI')
        hdu = fits.PrimaryHDU(data=result, header=tgt_header)
        hdu.writeto(repr480_image_filenames[filtername], overwrite=True)


# In[13]:


for filtername in image_sub_filenames:
    if not os.path.exists(repr480_image_sub_filenames[filtername]):
        print(f"Reprojecting {filtername} {image_sub_filenames[filtername]} to {repr480_image_sub_filenames[filtername]}")
        result,_ = reproject.reproject_interp(image_sub_filenames[filtername], tgt_header, hdu_in='SCI')
        hdu = fits.PrimaryHDU(data=result, header=tgt_header)
        hdu.writeto(repr480_image_sub_filenames[filtername], overwrite=True)


# In[14]:


rgb = np.array([
    fits.getdata(repr480_image_filenames['f210m']),
    fits.getdata(repr480_image_filenames['f300m']),
    fits.getdata(repr480_image_filenames['f360m'])
]).swapaxes(0,2).swapaxes(0,1)


# In[15]:


rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,0]),
                       simple_norm(rgb[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,1]),
                       simple_norm(rgb[:,:,2], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
plt.figure(figsize=(24,10))
plt.imshow(rgb_scaled, origin='lower')
plt.xticks([]);
plt.yticks([]);


# In[15]:


# rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,0]),
#                        simple_norm(rgb[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,1]),
#                        simple_norm(rgb[:,:,2], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
# plt.figure(figsize=(24,10))
# plt.imshow(rgb_scaled, origin='lower')
# plt.xticks([]);
# plt.yticks([]);


# In[16]:


png_path = '/orange/adamginsburg/jwst/sgrb2/pngs'
save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_210-300-360.png', avm=AVM)


# In[17]:


rgb2 = np.array([
    fits.getdata(repr480_image_filenames['f210m']),
    fits.getdata(repr480_image_filenames['f300m']),
    fits.getdata(repr480_image_filenames['f480m'])
]).swapaxes(0,2).swapaxes(0,1)


# In[18]:


rgb_scaled2 = np.array([simple_norm(rgb2[:,:,0], stretch='asinh', min_percent=1, max_percent=99)(rgb2[:,:,0]),
                       simple_norm(rgb2[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb2[:,:,1]),
                       simple_norm(rgb2[:,:,2], stretch='asinh', min_percent=1, max_percent=99)(rgb2[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
plt.figure(figsize=(24,10))
plt.imshow(rgb_scaled2, origin='lower')
plt.xticks([]);
plt.yticks([]);


# In[19]:


rgb_scaled.shape


# In[20]:


# save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_210-300-480.png', avm=AVM)


# In[21]:


# rgb3 = np.array([
#     fits.getdata(repr480_image_filenames['f480m']),
#     fits.getdata(repr480_image_filenames['f300m']),
#     fits.getdata(repr480_image_filenames['f210m']),
# ]).swapaxes(0,2).swapaxes(0,1)
# save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_480-300-210.png', avm=AVM)


# In[ ]:





# In[ ]:





# In[22]:


F150W = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f150w/jw05365-o001_t001_nircam_clear-f150w_i2d.fits')
F182M = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f182m/jw05365-o001_t001_nircam_clear-f182m_i2d.fits')
F187N = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f187n/jw05365-o001_t001_nircam_clear-f187n_i2d.fits')
F210M = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f210m/jw05365-o001_t001_nircam_clear-f210m_i2d.fits')
F212N = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f212n/jw05365-o001_t001_nircam_clear-f212n_i2d.fits')
F300M = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f300m/jw05365-o001_t001_nircam_clear-f300m_i2d.fits')
F360M = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f360m/jw05365-o001_t001_nircam_clear-f360m_i2d.fits')
F405N = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_f405n-f444w/jw05365-o001_t001_nircam_f405n-f444w_i2d.fits')
F410M = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f410m/jw05365-o001_t001_nircam_clear-f410m_i2d.fits')
F466N = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_f444w-f466n/jw05365-o001_t001_nircam_f444w-f466n_i2d.fits')
F480M = fits.open('/orange/adamginsburg/jwst/sgrb2/mastDownload/JWST/jw05365-o001_t001_nircam_clear-f480m/jw05365-o001_t001_nircam_clear-f480m_i2d.fits')


# In[23]:


import reproject


# In[ ]:





# In[ ]:





# In[24]:


Bra = fits.open('/orange/adamginsburg/jwst/sgrb2/NB/F405_minus_F410cont.fits')


# In[25]:


Bra['SCI'].data.shape


# In[26]:


F480M['SCI'].data.shape


# In[27]:


F410M['SCI'].data.shape


# In[28]:


F360M['SCI'].data.shape


# In[29]:


F480data = F480M['SCI'].data[:2507,:5649]
F410data = F410M['SCI'].data[:2507,:5649]
F300data = F300M['SCI'].data[:2507,:5649]
Bradata = Bra['SCI'].data[:2507,:5649]
F360data = F360M['SCI'].data[:2507,:5649]
F466data = F466N['SCI'].data[:2507,:5649]
#F405data = F405N['SCI'].data[:2507,:5650]


# In[30]:


rgb = np.array([F466data, F410data, F300data]).swapaxes(0,2).swapaxes(0,1)


# In[31]:


# rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=1, max_percent=98)(rgb[:,:,0]),
#                        simple_norm(rgb[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,1]),
#                        simple_norm(rgb[:,:,2], stretch='asinh', min_percent=1, max_percent=98)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
# plt.figure(figsize=(24,10))
# plt.imshow(rgb_scaled, origin='lower')
# plt.xticks([]);
# plt.yticks([]);


# In[32]:


#save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_466-410-300.png', avm=AVM)


# In[33]:


rgb = np.array([F480data, F410data, F360data]).swapaxes(0,2).swapaxes(0,1)


# In[ ]:





# In[34]:


# rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=1, max_percent=98.5)(rgb[:,:,0]),
#                        simple_norm(rgb[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,1]),
#                        simple_norm(rgb[:,:,2], stretch='asinh', min_percent=1, max_percent=98.5)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
# plt.figure(figsize=(12,5))
# plt.imshow(rgb_scaled, origin='lower')
# plt.xticks([]);
# plt.yticks([]);


# In[35]:


save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_480-410-360.png', avm=AVM)


# In[40]:


F480data.shape, fits.getdata(repr480_image_filenames['f187n']).shape


# In[58]:


rgb = np.array([fits.getdata(repr480_image_filenames['f480m']),
                fits.getdata(repr480_image_filenames['f405n']),
                fits.getdata(repr480_image_filenames['f187n'])]).swapaxes(0,2).swapaxes(0,1)
save_rgb(rgb/np.nanmedian(rgb), f'{png_path}/SgrB2_RGB_480-405-187.png', avm=AVM)


# In[62]:


rgb = np.array([fits.getdata(repr480_image_filenames['f480m']),
                fits.getdata(repr480_image_filenames['f405n']),
                fits.getdata(repr480_image_filenames['f187n'])]).swapaxes(0,2).swapaxes(0,1)
rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=1, max_percent=97)(rgb[:,:,0]),
                       simple_norm(rgb[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,1]),
                       simple_norm(rgb[:,:,2], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_480-405-187_scaled.png', avm=AVM)


# In[43]:


# list of filters from long to short so it's in RGB order
filternames = list(image_filenames.keys())[::-1]
for f1, f2, f3 in zip (filternames, filternames[1:], filternames[2:]):
    print(f1,f2,f3)
    rgb = np.array([
        fits.getdata(repr480_image_filenames[f1]),
        fits.getdata(repr480_image_filenames[f2]),
        fits.getdata(repr480_image_filenames[f3]),
    ]).swapaxes(0,2).swapaxes(0,1)
    rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,0]),
                           simple_norm(rgb[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,1]),
                           simple_norm(rgb[:,:,2], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)

    save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_{f1[1:4]}-{f2[1:4]}-{f3[1:4]}.png', avm=AVM)


# In[44]:


filternames_sub = list(image_sub_filenames.keys())[::-1]
for f1, f2, f3 in zip (filternames_sub, filternames_sub[1:], filternames_sub[2:]):
    print(f1,f2,f3)
    rgb = np.array([
        fits.getdata(repr480_image_sub_filenames[f1]),
        fits.getdata(repr480_image_sub_filenames[f2]),
        fits.getdata(repr480_image_sub_filenames[f3]),
    ]).swapaxes(0,2).swapaxes(0,1)
    rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=5, max_percent=99)(rgb[:,:,0]),
                           simple_norm(rgb[:,:,1], stretch='asinh', min_percent=5, max_percent=99)(rgb[:,:,1]),
                           simple_norm(rgb[:,:,2], stretch='asinh', min_percent=5, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)

    save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_{f1}_{f2}_{f3}.png', avm=AVM)


# In[45]:


rgb.shape, rgb_scaled.shape


# In[60]:


rgb = np.array([fits.getdata(repr480_image_filenames['f480m']),
                fits.getdata(repr480_image_sub_filenames['f405n-f410m']),
                fits.getdata(repr480_image_sub_filenames['f187n-f182m'])]).swapaxes(0,2).swapaxes(0,1)
save_rgb(rgb/np.nanmedian(rgb), f'{png_path}/SgrB2_RGB_480-405m410-187m182.png', avm=AVM)


# In[63]:


rgb = np.array([fits.getdata(repr480_image_filenames['f480m']),
                fits.getdata(repr480_image_sub_filenames['f405n-f410m']),
                fits.getdata(repr480_image_sub_filenames['f187n-f182m'])]).swapaxes(0,2).swapaxes(0,1)
rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=1, max_percent=97)(rgb[:,:,0]),
                       simple_norm(rgb[:,:,1], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,1]),
                       simple_norm(rgb[:,:,2], stretch='asinh', min_percent=1, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_480-405m410-187m182_scaled.png', avm=AVM)


# In[47]:


repr480_image_sub_filenames


rgb = np.array([fits.getdata(repr480_image_sub_filenames['f480m-f360m']),
                fits.getdata(repr480_image_sub_filenames['f405n-f410m']),
                fits.getdata(repr480_image_sub_filenames['f187n-f182m'])]).swapaxes(0,2).swapaxes(0,1)
rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='asinh', min_percent=5, max_percent=99)(rgb[:,:,0]),
                       simple_norm(rgb[:,:,1], stretch='asinh', min_percent=5, max_percent=99)(rgb[:,:,1]),
                       simple_norm(rgb[:,:,2], stretch='asinh', min_percent=5, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_480m360-405m410-187m182_scaled.png', avm=AVM)
rgb_scaled = np.array([simple_norm(rgb[:,:,0], stretch='log', min_percent=5, max_percent=99)(rgb[:,:,0]),
                       simple_norm(rgb[:,:,1], stretch='log', min_percent=5, max_percent=99)(rgb[:,:,1]),
                       simple_norm(rgb[:,:,2], stretch='log', min_percent=5, max_percent=99)(rgb[:,:,2])]).swapaxes(0,2).swapaxes(0,1)
save_rgb(rgb_scaled, f'{png_path}/SgrB2_RGB_480m360-405m410-187m182_logscaled.png', avm=AVM)
