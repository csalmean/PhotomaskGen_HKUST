# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 14:42:15 2021

Greyscale photolithography requires the breaking down of each shape into many small (sub-wavelength) pixels. Each pixel has a window size of which can be
used to adjust the amount of transmitted light.
By this means, we control the light penetration depth, allowing the development of 3D microstructures.

Since each pixel must be smaller than the resolution of the mask aligner, many millions of pixels must be defined. This can't be done graphically
so this code allows the photomasks to be constructed hierarchically using the following code.

This program works as follows:
    - Wafer conditioning: insertion of the alignment marks, protective borders and
     hard etch-masks. Alignment marks are awful as I am stuck with the heater design
     from a previous project. We also only use 21/27 chips as arraying the remaining 6 is too complicated.
     
    - Input equations and parameters for the desired 3D shape.
    - Generation of table representing every pixel in x/y, with height as the parameter held in each cell.
    - Change from height to pixel size.
    - Generate pixel array at this point.

@author: Chris Salmean
"""
import numpy as np
import pandas as pd
import math
import time
from scipy.ndimage import interpolation
from scipy import ndimage, misc
from PIL import Image
from collections import defaultdict
import gdspy

# Chips are designed for BRIGHT FIELD MASK. However, tests showed best results on dark-field shapes.
# We therefore include borders around each chip and around the entire chip array.
# This will mean small pixels = shallow etching, large pixels = deep etching

# timing calcualations:
t=time.time()

# First we determine the patterned area characteristics. 
# Dimensions of entire chip (including its border and the buffer gap) in microns
chipxdimension=6500
chipydimension=24000

# Patterned area dimensions
patterned_xdim=2000 #um
patterned_ydim=10000 #um

# Clear area on either side of patterned area
xbuffer=0 #um
ybuffer=6000 #um

# Space inbetween each chip (i.e. array spacing)
x_chip_spacing=500 #um
y_chip_spacing=100 #um

# Outline chip array dimensions on wafer
wafer_height=3 #chips
wafer_width=7 #chips

# Initialise DGSII library to contain shapes.
lib=gdspy.GdsLibrary()
gdspy.current_library = gdspy.GdsLibrary()

writer=gdspy.GdsWriter('FB greyscale Final.gds')
wafer=lib.new_cell('wafer')

# Mask Conditioning (squares, protector, alignment marks)
# Squares for mask writer calibration:
box=gdspy.Rectangle((0,0),(1,1))
boxcell=lib.new_cell('boxcell').add(box)
cond=lib.new_cell('cond')
cond.add(boxcell)
cond.add(boxcell.copy('box2',translation=(99999,99999)))

# Area of the wafer which will be protected from any patterning in BF configuration
protection=gdspy.Rectangle((10,10),(99990,99990))

# Now we want to make a clearing in the centre for the etched chips.
clearing_xcorner=50000-(wafer_width*0.5*(chipxdimension))
clearing_ycorner=47436.2-(wafer_height*0.5*(chipydimension))
clearing=gdspy.Rectangle((clearing_xcorner,clearing_ycorner),(100000-clearing_xcorner,94872.4-clearing_ycorner))

# Because of the (admittedly poorly-judged) placement of the heater alignment marks, backside alignment is
# difficult. Instead, we should define where exactly the heaters are meant to be.

# Heater outline for chips 1 and 13 on the middle row.
# We need to specify exactly where these boxes go:
clearing2=gdspy.Rectangle((10000,42318.2),(12000,52554.2))
clearing3=gdspy.Rectangle((88000,42318.2),(90000,52554.2))
combinedclearing=gdspy.boolean(clearing,clearing2,'or')
combinedclearing2=gdspy.boolean(combinedclearing,clearing3,'or')

outerborder=lib.new_cell('outerborder').add(gdspy.boolean(protection,combinedclearing2,'not'))

# Alignment marks
# Horrendous code but the marks on my previous masks were unneccessarily complex. In future iterations, it would be best to simplify this.
# The marks are mirrored twice. Simply give position of first mark (top left in this case)
xpos=22170.04
ypos=64516.16

alignmentmark=lib.new_cell('AM')

try:
    # Small crosses first:
    arm_length=102
    arm_thickness=4
    
    temp_xloc=xpos-150
    temp_yloc=ypos+97
    
    arm1=gdspy.Rectangle((temp_xloc,temp_yloc),(temp_xloc+arm_length,temp_yloc+arm_thickness))
    
    temp_xloc+=int(0.5*(arm_length-arm_thickness))
    temp_yloc-=int(0.5*(arm_length-(1.5*arm_thickness)))
    arm2=gdspy.Rectangle((temp_xloc,temp_yloc),(temp_xloc+arm_thickness,temp_yloc+arm_length))
    
    cross=gdspy.boolean(arm1,arm2,'or')
    alignmentmark.add(cross)
    
    # Now make bar arrays   
    length=232
    width=12
    number=12
    
    temp_yloc=ypos-0.96
    temp_xloc=xpos-481.04
    barloc=(temp_xloc,temp_yloc)
    
    bar=lib.new_cell('bar')
    bar.add(gdspy.Rectangle((0,0),(length,width)))
    spacing=(0,(2*(width-2)))
    bararray1=gdspy.CellArray(bar,1,number,spacing,barloc)      
    
    AM=[]
    AM.extend(bararray1.get_polygonsets())
    
    bararray=bararray1.get_polygonsets()
    
    # Make mirrored versions of these arrays.
    for poly in bararray:
        poly.mirror((xpos,ypos),((xpos-1),(ypos+1)))
        AM.append(poly)   
        
    # Now make squared crosses. Perhaps best just to take a series of points to
    # Define 1 arm then copy with rotation?
    # This should all be much simpler, but we are unfortunately prisoners of previous mask design errors.
    xsquare=gdspy.Rectangle((21689,64855.2),(21741,64907.2))
    xsquare2=gdspy.Rectangle((21739,64867.7),(21781,64894.7))
    xsquare3=gdspy.Rectangle((21780,64875.2),(21793.5,64887.2))
    xsquare4=gdspy.boolean(xsquare3,xsquare2,'or')
    xsquare5=gdspy.Rectangle((21780,64877.7),(21805,64884.7))
    xsquare6=gdspy.boolean(xsquare5,xsquare4,'or')
    squarecrossarm=gdspy.boolean(xsquare,xsquare6,'or')
    squarecrossarm2=gdspy.copy(squarecrossarm).mirror((21805,64881.2),(21805,64880.2))
    squarecrossarm3=gdspy.copy(squarecrossarm2).rotate((np.pi/2),(21805,64881.2))
    squarecrossarm4=gdspy.copy(squarecrossarm3).mirror((21805,64881.2),(21806,64881.2))
    
    squarecrosspc=gdspy.boolean(squarecrossarm,squarecrossarm2,'or')
    squarecrosspc2=gdspy.boolean(squarecrossarm3,squarecrossarm4,'or')

    AM.append(gdspy.boolean(squarecrosspc,squarecrosspc2,'or'))
    
    # Make scaled versions of the bars/cross.
    for i in range (len(AM)):
        AM.append(gdspy.copy(AM[i]).scale(2,2,(xpos,ypos)))
        AM.append(gdspy.copy(AM[i]).scale(4,4,(xpos,ypos)))
   
    AM.append(cross) 
    
    # Now we have a full alignment mark, mirror twice
    
    for i in range (len(AM)):
        AM.append(gdspy.copy(AM[i]).mirror((50000,0),(50000,1)))
    for i in range (len(AM)):
        AM.append(gdspy.copy(AM[i]).mirror((0,47436.21),(1,47436.21)))     
    
    for i in range (len(AM)):
        alignmentmark.add(AM[i])
  
except:
    print('Failed to make alignment marks') 

# Now just boolean to cut these marks out of the border
cond.add(gdspy.boolean(outerborder,alignmentmark,'not'))
wafer.add(cond)

# Next we define the borders of each chip.
# Best to follow dimensions of old chips; 2x r1 circles with 20mm centre-centre
outer_border_xdimension=chipxdimension-(x_chip_spacing)
outer_border_ydimension=chipydimension-(y_chip_spacing)

#define cutout shape
endpc_radius=(0.5*patterned_xdim)+xbuffer
separating_distance=18000

# Generate list of chip locations (beginning in southwest corner) and place borders
borders=lib.new_cell('borders')
locations=[]
yloc=clearing_ycorner
xloc=clearing_xcorner

for j in range(wafer_width):
    locations.append((xloc,yloc))
    
    outer_border=gdspy.Rectangle((xloc,yloc),(xloc+outer_border_xdimension,yloc+outer_border_ydimension))
    outer_border.translate((0.5*x_chip_spacing), (0.5*y_chip_spacing))
    
    inner_border_e1=gdspy.Round((xloc,yloc),endpc_radius)
    inner_border_e2=gdspy.Round((xloc,yloc+separating_distance),endpc_radius)
    inner_border_ec=gdspy.boolean(inner_border_e1,inner_border_e2,'or')
    inner_border_c=gdspy.Rectangle((xloc-(xbuffer+(0.5*(patterned_xdim))),yloc),
                                   (xloc+(xbuffer+(0.5*(patterned_xdim))),(separating_distance+yloc)))

    inner_border_combined=gdspy.boolean(inner_border_ec,inner_border_c,'or') 
    inner_border_combined.translate((chipxdimension/2),3000)
    
    borders.add(gdspy.boolean(outer_border,inner_border_combined,'not'))
    
    for i in range((wafer_height-1)):
        yloc+=chipydimension
        locations.append((xloc,yloc))
        outer_border=gdspy.Rectangle((xloc,yloc),(xloc+outer_border_xdimension,yloc+outer_border_ydimension))
        outer_border.translate((0.5*x_chip_spacing), (0.5*y_chip_spacing))
        inner_border_e1=gdspy.Round((xloc,yloc),endpc_radius)
        inner_border_e2=gdspy.Round((xloc,yloc+separating_distance),endpc_radius)
        inner_border_ec=gdspy.boolean(inner_border_e1,inner_border_e2,'or')
        inner_border_c=gdspy.Rectangle((xloc-(xbuffer+(0.5*(patterned_xdim))),yloc),
                                       (xloc+(xbuffer+(0.5*(patterned_xdim))),(separating_distance+yloc)))
    
        inner_border_combined=gdspy.boolean(inner_border_ec,inner_border_c,'or') 
        inner_border_combined.translate((chipxdimension/2),3000)
        
        borders.add(gdspy.boolean(outer_border,inner_border_combined,'not'))
    
    yloc=clearing_ycorner
    xloc+=chipxdimension

# wafer.add(borders)

tt=time.time()
print('Wafer conditioned in '+str(tt-t))

# Next for the important part: generation and printing of greyscale patterns. Starting from unit cell bitmaps, we form an array with 
# correct width and height to furnish a single chip. We then place in arrays using the list of chip locations as each chip is copied three times.

# Create list of parameters to feed into generators
pixelpitch=3.2 #grey levels: 12
minimum_pixel_width=1.0 #um
machine_grid=0.2 #um

# Since we are designing a range of shapes, some are best defined in cartesian
# terms while others are better in polar terms. We need to treat them separately.

# First create a unit cell's pixel map (using dimensions and pixel size)
def unitmapgenerator(x_size,y_size,pixel_pitch):
    # create structure unit cell then stitch together afterwards
    x_size/=pixel_pitch
    y_size/=pixel_pitch
    
    # Create array with same number of rows and columns as there are pixels
    # in the unit cell.
    
    unit_map=np.ndarray((int(x_size),int(y_size)))
    unit_map*=0
    unit_map=abs(np.nan_to_num(unit_map,nan=0))
    return unit_map

# For ramps and pleateaux we need to determine proportion of distance to line where h=1.
# We then use this to determine height at this point. Only ened to do 0deg and 90deg

def cartesian_builder(unit_map,pixelpitch,rotation_angle,after_rot,shape,plateau_prop,gap_prop,direction):      
    # IF current location on slope (sqrt(x**2 + y**2)) <=slope length (= unit length - plateau - blank)
    # We can use h=ml = sqrt(i**2 + j**2)/slope_length
    # ELIF location on slope within plateau area, h=1,
    # ELSE h=0
    
    # This sqrt(x**2 + y**2) actually comes from trig; isin(dir)+jcos(dir)
    rotation_angle_rad=rotation_angle*np.pi/180
    
    if shape=="ramp" or shape=="plateau" or shape=="thorn":
        total_vector_length=np.size(unit_map,0)*math.cos(rotation_angle_rad)+np.size(
            unit_map,1)*abs(math.sin(rotation_angle_rad))
        slope_length=(1-plateau_prop-gap_prop)*total_vector_length
        
        for i in range(np.size(unit_map,0)):
            for j in range(np.size(unit_map,1)): 
                current_length=i*math.cos(rotation_angle_rad)+j*abs(math.sin(rotation_angle_rad))
                
                # height is equal to length * slope 
                # apply to determine HEIGHT of shape at each location   
                unit_map[i,j]= (current_length)/(slope_length)
                if (current_length>=slope_length)and(current_length<=(total_vector_length*(1-gap_prop))):
                    unit_map [i,j]=1
                elif(current_length>slope_length)&(current_length>(total_vector_length*(1-gap_prop))):
                    unit_map[i,j]=0

    elif shape=="NACA":
        #  r = sqrt(((1+y)(1-y)^3)/8)
        #  h = sqrt (r**2 - i**2)
        # NACA 00xx:
            # r=5[xx]*(0.2969y^0.5 - 0.1260y - 0.3516y^2 + 0.2843y^3 - 0.1015y^4)
        
        # Isolate midpoint of the grid first:
        midx=np.size(unit_map,0)/2
        # midy=np.size(unit_map,1)/2
        midy=0
        shape_length=np.size(unit_map,1)*(1-gap_prop)/np.size(unit_map,1)
        NACA_thickness = 0.5
        
        for i in range(np.size(unit_map,0)):
            for j in range(np.size(unit_map,1)):
                
                proportionalx=(i-midx)/(np.size(unit_map,0))/shape_length
                proportionaly=(j-midy)/(np.size(unit_map,1))/shape_length
                
                # Apply equation from above to determine HEIGHT of shape at each location
                r = 5*NACA_thickness*(0.2969*(proportionaly**0.5) - 
                            0.1260*proportionaly -
                            0.3516*proportionaly**2 + 
                            0.2843*proportionaly**3 - 
                            0.1015* proportionaly**4)
                
                # r = np.sqrt((1+proportionaly)*((1-proportionaly)**3)/8)
                if abs(proportionalx)>r:
                    unit_map[i,j]=0
                else: 
                    unit_map[i,j]=np.sqrt(r**2 - proportionalx**2)     
                    
        unit_map[np.isnan(unit_map)]=0
       
    if after_rot>0:
        unit_map=ndimage.rotate(unit_map, after_rot, reshape=True,order=1)
        unit_map/=np.max(abs(unit_map))
    
    if direction=="indent":
        unit_map=1-unit_map  
    
    return unit_map

# For polar we want to indicate central axis of shape and define from there.
def polar_builder(unit_map,pixelpitch,shape,R,torus_Rm,direction):
    midx= np.size(unit_map,0)/2
    midy= np.size(unit_map,1)/2              
    
    # Sphere, hemitorus
    # sphere with centre (mpx,mpy), radius R, projected vector length v.
    # surface is h=sqrt(R^2-v^2)
    # h = R^2 - (i-mpx)^2 - (j-mpy)^2
    
    # torus with centre (mpx, mpy), major radius R and minor radius Rm
    # equation is h=sqrt(Rm**2 - v**2)
    # vector v is relative to midpoint of minor circle.
    # v = sqrt((i-mpx)**2 + (j-mpy)**2)-(R)

    R/=pixelpitch
    torus_Rm/=pixelpitch
    
    for i in range(np.size(unit_map,0)):
        for j in range(np.size(unit_map,1)):
            ilen=abs(midx-i)
            jlen=abs(midy-j)
            v=np.sqrt(ilen**2+jlen**2)
            
            if shape=="sphere":
                if v>R:
                    pass
                else:
                    #rs^2 = rc^2 + (z-h)^2
                    unit_map[i,j]=(R**2 - v**2)**0.5
            
            elif shape=="torus":
                
                v=np.sqrt((i-midx)**2+(j-midy)**2)-(R)
                
                if (abs(v)>torus_Rm):
                    pass
                else:      
                    # rm^2 = (rc-rs-i)^2 + (rc-rs-j)^2 +(z-h)^2
                    unit_map[i,j]=(R**2-v**2)
                    
            
    unit_map/=np.max(abs(unit_map))
    
    if direction=="indent":
        unit_map=1-unit_map    
    return unit_map

# To form greyscale image for pre-inspection
# gsmap=unit_map.copy()
# gsmap*=255
# gsmap=gsmap.round(0)
# im=Image.fromarray(gsmap)
# im = im.convert("L")
# im.save('test3.png')

designs=defaultdict(list)

x_size=200 #um
y_size=200 #um

#specify coordinate system first with "c" or "p"
# ramp with 150um length and 100um width, rotated to 0 degrees

# Ramps
designs[0]=({
    "grid":"c",
    "rot":90,
    "after_rot":180,
    "shape":"ramp",
    "plateau_prop":0,
    "gap_prop":0,
    "dir":"raise"})

# Plateaux
designs[1]=({
    "grid":"c",
    "rot":90,
    "after_rot":180,
    "shape":"plateau",
    "plateau_prop":0.25,
    "gap_prop":0.25,
    "dir":"raise"})

# Dimples
designs[2]=({
    "grid":"p",
    "after_rot":0,
    "shape":"sphere",
    "R":50,
    "Rm":0,
    "dir":"indent"})

# Ring cavities x2
designs[3]=({
    "grid":"p",
    "shape":"torus",
    "after_rot":0,
    "R":35,
    "Rm":15,
    "dir":"indent"})

designs[4]=({
    "grid":"p",
    "shape":"torus",
    "after_rot":0,
    "R":25,
    "Rm":25,
    "dir":"indent"})

# Teardrops
designs[5]=({
    "grid":"c",
    "rot":0,
    "after_rot":180,
    "shape":"NACA",
    "NACA":"0050",
    "plateau_prop":0,
    "gap_prop":0.4,
    "dir":"raise"})


# Thorns
designs[6]=({
    "grid":"c",
    "rot":225,
    "after_rot":315,
    "shape":"thorn",
    "plateau_prop":0,
    "gap_prop":0.5,
    "dir":"raise"})
        
# We generate pixel grid for a single unit cell.
for variation,design in (designs.items()):
    designs[variation]['unitmap']=[]
    print("Design: "+str(variation))
    
    # We generate pixel grid for a single unit cell.
    unit_map=unitmapgenerator(x_size,y_size,pixelpitch)    
    
    # Treat differently, depending if cartesian or polar.
    if "c" in designs[variation]["grid"]:
        
        unit_map=cartesian_builder(unit_map,pixelpitch,
                                   designs[variation]["rot"],
                                   designs[variation]["after_rot"],
                                   designs[variation]["shape"],
                                   designs[variation]["plateau_prop"],
                                   designs[variation]["gap_prop"],
                                   designs[variation]["dir"])
        
    elif "p" in designs[variation]["grid"]:
        
        unit_map=polar_builder(unit_map,pixelpitch,
                               designs[variation]["shape"],
                               designs[variation]["R"],
                               designs[variation]["Rm"],
                               designs[variation]["dir"])
                   
    unit_map/=np.max(unit_map)
    designs[variation]['unitmap']=unit_map
    #Take a copy and normalise to 255 bits for production of greyscale image
    gsmap=unit_map.copy()
    gsmap*=255
    gsmap=gsmap.round(0)
    gsmap=255-gsmap
    
    im=Image.fromarray(gsmap)
    im = im.convert("L")
    im.save(designs[variation]["shape"]+'.png')
    
    

ttt=time.time()
print('Generated unit pixmaps in '+str(ttt-tt))

# Now generate the pixmaps with pixel size depending upon height in unit map
def gdsii_array_builder(pixmap,pitch,chipnumber,field, after_rot,direction):
    
    # First we must determine if the pixels should be larger for deep etch or smaller.
    # If BF, means pixel array is open and pixels are larger where eching should be less shallow.
    # In DF, pixel array is covered and pixels give WINDOWS so are smaller where etching should be less shallow. 
    
    if field=='BF':
        pixel_range=pitch-(minimum_pixel_width)
    elif field=='DF':
        pixel_range=pitch-minimum_pixel_width+0.2
        pixmap=(1-pixmap)
    
    # First we determine the pixel fill dimensions:

    pixmap*=5*pixel_range
    pixmap=pixmap.round(0)
    pixmap/=5
    
    pixmap+=(minimum_pixel_width-0.2)
    pixmap[pixmap<minimum_pixel_width]=0
    
    if field=='DF':
        pixmap[pixmap>pitch]=pitch
    # um_map*=(pitch/np.max(um_map))

    gdsii_unit=gdspy.Cell('chip{} unit'.format(chipnumber),True)
    
    for i in range (np.size(pixmap,0)):
        for j in range(np.size(pixmap,1)):
            size=(pixmap[i,j],pixmap[i,j])
            origin=((i*pitch)+(0.5*(pitch-size[0])),(j*pitch)+(0.5*(pitch-size[1])))
            
            if (size[0]<minimum_pixel_width)or (size[1]<minimum_pixel_width):
                pass
            else:
                pixel=gdspy.Rectangle(origin,(origin[0]+size[0],origin[1]+size[1]))
                gdsii_unit.add(pixel)


    # Now that we have a unit cell, we invert if the cell is for DF
    final_unit=gdspy.Cell('completed_unitcell_{}'.format(chipnumber),True)

    if field=='DF':
        BG_xdim=(np.size(pixmap,0)*pitch)
        BG_ydim=(np.size(pixmap,1)*pitch)
        BGfill=gdspy.Rectangle((0,0),(BG_xdim,BG_ydim))
        
        final_unit.add(gdspy.boolean(BGfill,gdsii_unit,'not'))
    
    elif field=='BF':
        final_unit=gdsii_unit
    
    # wafer.add(final_unit)
    # We must avoid pixels being placed directly next to each other by the array.
    
    gdsii_array=gdspy.Cell('chip{} array'.format(chipnumber),True)
    
    #since we are multiplying up by a decimal the scaled pixmap may be a couple of microns short.
    # we need to calculate the error so the array can be perfectly centred.
    after_rot*=np.pi/180
        
    spacing=((np.size(pixmap,0)*pitch),(pitch*np.size(pixmap,1)))
    arraycolumns=int(patterned_xdim/spacing[0])
    arrayrows=int(patterned_ydim/spacing[1])

    # we want to place each array in triplicate
    for i in range(3):
        
        infill=gdspy.Cell('chip{} infill'.format(chipnumber),True)
        
        origin_x=locations[chipnumber][0]+(0.5*(chipxdimension-patterned_xdim))+8
        origin_y=locations[chipnumber][1]+(0.5*(chipydimension-patterned_ydim))+8
        origin=(origin_x,origin_y)
        
        gdsii_array.add(gdspy.CellArray(final_unit, arraycolumns, arrayrows, spacing,origin))
        
        # If the shapes are indented into a flat surface we also fill the channel before/after the heated area
        if direction=="indent":
            
            indentsize=(patterned_xdim,(ybuffer-(2*endpc_radius)))
            origin=((origin_x-8),(origin_y-indentsize[1]))
            gdsii_array.add(gdspy.Rectangle(origin,(origin[0]+indentsize[0],origin[1]+indentsize[1])))            

            origin=((origin_x-8),(origin_y+(spacing[1]*arrayrows)))
            gdsii_array.add(gdspy.Rectangle(origin,(origin[0]+indentsize[0],origin[1]+indentsize[1])))   
        
        chipnumber += 1
    
    writer.write_cell(final_unit)
    wafer.add(infill)
    wafer.add(gdsii_array)
    return pixmap

def textgenerator(layout,chipnumber):
    if 'NACA' in layout.keys():
        layout['shape']+=layout['NACA']
        
    for i in range(3):
           
        text=gdspy.Cell('chip{} text'.format(chipnumber),True)
        origin=(locations[chipnumber][0]+500,locations[chipnumber][1]+500)
        
        if layout['grid']=='c':
            string=(layout['shape']+' '+str(chipnumber)+', p:'+str(layout['plateau_prop'])+', g:'+
                                              str(layout['gap_prop']))
        else:
            string=(layout['shape']+' '+str(chipnumber)+', R:'+str(layout['R'])+', Rm:'+
                                              str(layout['Rm']))
        
        label=gdspy.Text(string,225,origin)
        mirrorx=(locations[chipnumber][0]+(0.5*chipxdimension))
        mirrory=(locations[chipnumber][1]+(0.5*chipydimension))
        
        text.add(gdspy.PolygonSet.mirror(label,(mirrorx,mirrory),(mirrorx,mirrory+1)))
        
        alltext.add(text)
        chipnumber+=1
        
# Now we want to print each design in triplicate, moving across the wafer.

alltext=gdspy.Cell('Text')

for design,layout in enumerate(designs.values()):
    print('Row:'+str(design))
    
    chipnumber=(design*wafer_height)
    chipcell=gdspy.Cell('Chip{}'.format(chipnumber))

    gdsii_array_builder(layout['unitmap'], pixelpitch,chipnumber,'DF',layout['after_rot'],layout['dir'])
 
    textgenerator(layout,chipnumber)
    

wafer.add(gdspy.boolean(borders, alltext, 'not'))      
tttt=time.time()
print('Generated gdsii arrays in '+str(tttt-ttt))
print('Done.')

# write the file now.
writer.write_cell(wafer)
writer.close()