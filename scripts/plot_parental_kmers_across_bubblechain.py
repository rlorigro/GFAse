import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib
import numpy as np
import scipy.stats as stats
import argparse


def fill_lists(input_file):
  num_kmers_in_path=0
  num_kmers_in_mat=0
  num_kmers_in_pat=0
  # store the kmer locations in x and y_pat, y_mat lists
  inFile=open(input_file,'r')

  # Read through the input file
  for line in inFile:
      splitLine=line.strip().split(',')

      # skip header line
      if splitLine[0] == "path_index":
        continue

      x_path_index_list.append(float(splitLine[0]))
      y_paternal_list.append(float(splitLine[1]))
      y_maternal_list.append(float(splitLine[2]))
      # if splitLine[1]=="0" and splitLine[2]=="0":
      #   x_no_match_list.append(int(splitLine[0]))
      #   y_no_match_list.append(int(splitLine[1]))
      
      # pat
      if splitLine[1]!="0":
        num_kmers_in_pat+=1
        x_pat_path_index_list.append(float(splitLine[0]))
        y_paternal_scatter_list.append(float(splitLine[1])) 
      # mat
      if splitLine[2]!="0":
        num_kmers_in_mat+=1
        x_mat_path_index_list.append(float(splitLine[0]))
        y_maternal_scatter_list.append(float(splitLine[2]))

      num_kmers_in_path+=1
  return num_kmers_in_path, num_kmers_in_mat, num_kmers_in_pat

def plot_data(panel, in_file,scale):
  # Plotting
  swarmY(x,center,ptwidth)
  y,x = swarmY(x_pat_path_index_list,scale,20)
  # print('x',x,'\ny',y)
  panela = panel.twiny()

  # plot paternal kmer locations
  panela.plot(x,y,#x_pat_path_index_list,y_pat_scatter_scaled,
              marker='o',
              markerfacecolor=(51/255,153/255,255/255),
              markeredgecolor='black',
              markersize=2,
              markeredgewidth=0,
              linewidth=0,
              )

  # plot maternal kmer locations
  panela.plot(x_mat_path_index_list,y_mat_scatter_scaled,
              marker='o',
              markerfacecolor=(243/255,35/255,209/255),
              markeredgecolor='black',
              markersize=2,
              markeredgewidth=0,
              linewidth=0,
              )

  panela.tick_params(bottom=True, labelbottom=True,
                     left=False, labelleft=False,
                     right=False, labelright=False,
                     top=False, labeltop=False)

  panel.set_xticks([int(x_path_index_list[0]),int(x_path_index_list[-1])])
  panel.set_xticklabels([int(x_path_index_list[0]),int(x_path_index_list[-1])])

  panel.set_xlabel("component kmers: " + str(num_kmers_in_path) + "\npat: "+ str(num_kmers_in_pat) + "   mat: " + str(num_kmers_in_mat),fontsize = 7 )

  panel.set_title("component: "+str(in_file.split("/")[-1][:-4]))

def windowed_sum_dif(panel,in_file,color):
  y_pat = np.convolve(np.array(y_paternal_list),np.ones(10000,dtype=int),'valid')

  y_mat = np.convolve(np.array(y_maternal_list),np.ones(10000,dtype=int),'valid')

  y_range = [abs(np.min(y_pat)),np.max(y_pat),abs(np.min(y_mat)),np.max(y_mat)]

  plot_max = np.max(y_range)
  plot_min = np.min(y_range)
  # y_plotting_list=y_pat - y_mat
  # print('ypat', y_pat, '\nymat', y_mat, '\ndif:',y_plotting_list)
  x = np.arange(0,len(y_pat))

  mat_col = (243/255,35/255,209/255)
  pat_col = (51/255,153/255,255/255)

  # color = [pat_col if y>0 else mat_col for y in y_plotting_list]

  # panela = panel.twinx()
  # panela = panel.twiny()
  panel.fill_between(x,y_pat,
                color=pat_col,
                linewidth=1,
                alpha=0.5
                )
  y_mat_neg = [y*-1 for y in y_mat]

  panel.fill_between(x,y_mat_neg,
                color=mat_col,
                linewidth=1,
                alpha=0.5
                )

  # plot 0 line
  panel.axhline(linewidth=1, color=(192/255,192/255,192/255))

  panel.set_title("component: "+str(in_file.split("/")[-1][:-4]), fontsize = 2)

  return float(plot_max),float(plot_min)

def distancePythag(x,xlist,ylistSqared):
  '''calculates distance using pythagorean theorem
         d = sqrt(xdist^2 + ydist^2)              '''
  lxl2 = np.power(xlist,2)           # square
  ladded = np.add(lxl2,ylistSqared)  # add
  ld = np.sqrt(ladded)               # square root
  lptdistance = [min(ld),x]         # min distnace to point
  return lptdistance

def swarmY(parent_kmers,center,ptwidth):
  ''' returns three lists x,y,qv of x asjusted points 
      from a SORTED list of two elements
      that fit within the ptwidth on either side

      ptsize=0.7                # point size in points
      ptDiameter = (ptsize)/72  # diameter of a point using .plot() in inches
      shift=ptDiameter/3        # set x shift as 1/3 of a point diameter
      x,y,c = swarmX(scs_q2,cov,0.45)
      '''
  # x,y,qv(color) lists
  xs=[]
  ys=[]
  # c=[]
  # plot the first points where they are
  xs.append(center)
  ys.append(parent_kmers[0])
  # c.append(sampId[0][1])
  # loop over the rest of the points moving them if they are 
  # close to other poings
  print('len(parent_kmers)',len(parent_kmers))
  for i in range(1,len(parent_kmers)):
    xcoor=center
    ycoor=parent_kmers[i]
    # for items less than 1 point diameter away adjust the x value
    if (((ycoor - ys[len(ys)-1])/25)*mainPanelHeight) < (ptDiameter):
      # print('this',(((ycoor - ys[len(ys)-1])/25)*mainPanelHeight),'ptDiameter',ptDiameter)

      placed=False
      lx = xcoor 
      rx = xcoor 
      # move the point to the left or right
      while placed==False:
        # keep track of pt distances on left and right
        rxl=[] # right x list
        lxl=[] # right x list
        y=[]   # list of y points
        j=len(xs)-1
        # only check the distance of points within a diameter on the y axis
        while ( ((ycoor - ys[j])/25)*mainPanelHeight) < (ptDiameter) and j>=0: 
          # print('move',lxl,((lx - xs[j])/12)*mainPanelWidth,rxl)
          lxl.append(((lx - xs[j])/12)*mainPanelWidth)
          rxl.append(((rx - xs[j])/12)*mainPanelWidth)
          y.append( ((ycoor-ys[j])/25)*mainPanelHeight)
          j-=1
        
        # distance calculations of nearby points
        y2 = np.power(y,2)
        # calc left point distances and increment 
        lptdistance = distancePythag(lx,lxl,y2)
        lx = lx - (shift)   # increment x position for next loop
        # calc right point distances and increment 
        rptdistance = distancePythag(rx,rxl,y2)#[min(rd),rx]
        rx = rx + (shift)  # increment x position for next loop

        # if either the left or right pt dist is > the diameter plot it there
        for dist in [rptdistance,lptdistance]:
          if dist[0]>ptDiameter:
            xs.append(dist[1])
            ys.append(parent_kmers[i])
            # c.append(sampId[i][1])
            placed=True
            break

        # drop points outside of x bounds
        if abs(lx - center) > ptwidth:
          placed=True
          break
        if abs(rx-center) > ptwidth:
          placed=True
          break
    else: 
      # plot isolated poitns in the center of the bin
      # if no points are nearby
        xs.append(center)
        ys.append(parent_kmers[i])
        # c.append(sampId[i][1])
  # switch order of these to make it swarm horizontal instead of vertical
  return xs,ys#,c

### parse args
parser=argparse.ArgumentParser()

parser.add_argument('-i','--input_file')
# parser.add_argument('-o','--output_file')

args=parser.parse_args()

input_file_1=args.input_file

# get the file name for the other haplotype path csv
file_ending = input_file_1[-4:]
other_hap = "0"
if input_file_1[:-4][-1]=="0":
  other_hap="1"

input_file_2 = input_file_1[:-5]+other_hap+file_ending

out_file = input_file_1[:-4]+"_"+other_hap+".png"

# set the figure size 
figureHeight=6
figureWidth=8

plt.figure(figsize=(figureWidth,figureHeight))


# set the size of the panels
  # Main Panels: 2'' wide and 2'' high
mainPanelHeight=2
mainPanelWidth=6

# set the variables for relative height and width of all 6 panels
relativeMainPanelWidth=mainPanelWidth/figureWidth
relativeMainPanelHeight=mainPanelHeight/figureHeight


### Main Panels
panel1= plt.axes([0.12,0.55,relativeMainPanelWidth,relativeMainPanelHeight])
panel2= plt.axes([0.12,0.08,relativeMainPanelWidth,relativeMainPanelHeight])

ptsize=2000                # point size in points
ptDiameter = (ptsize)/72  # diameter of a point using .plot() in inches
shift=ptDiameter/3        # set x shift as 1/3 of a point diameter


## Get data
# open the data file and store the log converted data into x and y lists

x_path_index_list=[]#np.array([])
y_paternal_list=[]#np.array([])
y_maternal_list=[]#np.array([])

y_paternal_scatter_list=[]
y_maternal_scatter_list=[]
x_pat_path_index_list=[]
x_mat_path_index_list=[]

##### Read through the 1st input file ##########
print('input_file_1', input_file_1)
num_kmers_in_path, num_kmers_in_mat, num_kmers_in_pat = fill_lists(input_file_1)

## Plotting
rel_max, rel_min = windowed_sum_dif(panel1,input_file_1,(243/255,35/255,209/255))
print('rel_max',rel_max,'rel_min',rel_min)

scale_plot =1
if abs(rel_min) > abs(rel_max):
  scale_plot = abs(rel_min)
  panel1.set_ylim(rel_min+(rel_min/3),0+abs(rel_min)+(abs(rel_min)/3))
  # panel1.set_yticks([rel_min+(rel_min/3),0+abs(rel_min)+(abs(rel_min)/3)])
  # since this gets changed - just have a # of kmers that aren't 0 and do [scale_plot]*kmers
  y_pat_scatter_scaled = [scale_plot for y in y_paternal_scatter_list ]
  y_mat_scatter_scaled = [scale_plot*-1 for y in y_maternal_scatter_list]
else:
  scale_plot = abs(rel_max)
  panel1.set_ylim(0-abs(rel_max)-(rel_max/3),rel_max+(rel_max/3))
  # panel1.set_yticks([0-abs(rel_max)-(rel_max/3),rel_max+(rel_max/3)])
  y_pat_scatter_scaled = [scale_plot for y in y_paternal_scatter_list ]
  y_mat_scatter_scaled = [scale_plot*-1 for y in y_maternal_scatter_list]

panel1.set_xlabel("component kmers: " + str(num_kmers_in_path) + "\npat: "+ str(num_kmers_in_pat) + "   mat: " + str(num_kmers_in_mat),fontsize = 7 )
panel1.set_title("component: "+str(input_file_1.split("/")[-1][:-4]))

x_path_index_list=[]
y_paternal_list=[]
y_maternal_list=[]

y_paternal_scatter_list=[]
y_maternal_scatter_list=[]
x_pat_path_index_list=[]
x_mat_path_index_list=[]

##### Read through the 2nd input file
print('\ninput_file_2', input_file_2)
num_kmers_in_path, num_kmers_in_mat, num_kmers_in_pat = fill_lists(input_file_2)

## Plotting
rel_max, rel_min = windowed_sum_dif(panel2,input_file_2,(0,0,0))
print('rel_max',rel_max,'rel_min',rel_min)
# y_pat_scatter_scaled = [rel_max *y for y in y_paternal_scatter_list ]
# y_mat_scatter_scaled = [rel_min *y for y in y_maternal_scatter_list]

# make plot balanced
if abs(rel_min) > abs(rel_max):
  scale = abs(rel_min)
  panel2.set_ylim(rel_min+(rel_min/3),0+abs(rel_min)+(abs(rel_min)/3))
  # panel2.set_yticks([rel_min+(rel_min/3),0+abs(rel_min)+(abs(rel_min)/3)])
  
  y_pat_scatter_scaled = [scale*y for y in y_paternal_scatter_list ]
  y_mat_scatter_scaled = [scale*y*-1 for y in y_maternal_scatter_list]
else:
  scale = abs(rel_max)
  panel2.set_ylim(0-abs(rel_max)-(rel_max/3),rel_max+(rel_max/3))
  # panel2.set_yticks([0-abs(rel_max)-(rel_max/3),rel_max+(rel_max/3)])
  y_pat_scatter_scaled = [scale*y for y in y_paternal_scatter_list ]
  y_mat_scatter_scaled = [scale*y*-1 for y in y_maternal_scatter_list]

# panel2.set_yticklabels(['mat','pat'])
panel2.set_xlabel("component kmers: " + str(num_kmers_in_path) + "\npat: "+ str(num_kmers_in_pat) + "   mat: " + str(num_kmers_in_mat),fontsize = 7 )
panel2.set_title("component: "+str(input_file_2.split("/")[-1][:-4]))


### Main Panel
panel1.tick_params(bottom=True, labelbottom=True,
                     left=True, labelleft=True,
                     right=False, labelright=False,
                     top=False, labeltop=False)
panel2.tick_params(bottom=True, labelbottom=True,
                     left=True, labelleft=True,
                     right=False, labelright=False,
                     top=False, labeltop=False)

plt.savefig(out_file,dpi=600)
plt.show()
