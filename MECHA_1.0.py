# -*- coding: utf-8 -*-

#Directory
dir='C:/Users/admintmp/Desktop/MECHA/'

#Inputs
Geom='PearlMillet2_Geometry.xml' #'Maize1_Geometry.xml'
Hydr='Hydraulics.xml'
BC='BC.xml'

#Libraries
import numpy as np #NumPy is the fundamental package for scientific computing with Python.
                   #It contains among other things:
                   #- a powerful N-dimensional array object
                   #- sophisticated (broadcasting) functions
                   #- tools for integrating C/C++ and Fortran code
                   #- useful linear algebra, Fourier transform, and random number capabilities
from numpy import genfromtxt #Load data from a text file, with missing values handled as specified.
from numpy.random import *  # for random sampling

import scipy.linalg as slin #Linear algebra functions

import pylab #Plotting functions library
from pylab import *  # for plotting

import networkx as nx #NetworkX is a Python language software package for the creation,
                      #manipulation, and study of the structure, dynamics, and functions of complex networks.
from lxml import etree #Tree element analysis module

import sys, os # On importe le module os qui dispose de variables 
               # et de fonctions utiles pour dialoguer avec votre 
               # système d'exploitation

#Import General data
OS=etree.parse(dir + 'in/General.xml').getroot().xpath('OS')[0].get("value")
Output_path=etree.parse(dir + 'in/General.xml').getroot().xpath('Output')[0].get("path")
Paraview=int(etree.parse(dir + 'in/General.xml').getroot().xpath('Paraview')[0].get("value"))
ParaviewWF=int(etree.parse(dir + 'in/General.xml').getroot().xpath('Paraview')[0].get("WallFlux"))
ParaviewMF=int(etree.parse(dir + 'in/General.xml').getroot().xpath('Paraview')[0].get("MembraneFlux"))
ParaviewPF=int(etree.parse(dir + 'in/General.xml').getroot().xpath('Paraview')[0].get("PlasmodesmataFlux"))
ParaviewWP=int(etree.parse(dir + 'in/General.xml').getroot().xpath('Paraview')[0].get("WallPot"))
ParaviewCP=int(etree.parse(dir + 'in/General.xml').getroot().xpath('Paraview')[0].get("CellPot"))
color_threshold=float(etree.parse(dir + 'in/General.xml').getroot().xpath('color_threshold')[0].get("value"))
thickness_disp=float(etree.parse(dir + 'in/General.xml').getroot().xpath('thickness_disp')[0].get("value"))
thicknessJunction_disp=float(etree.parse(dir + 'in/General.xml').getroot().xpath('thicknessJunction_disp')[0].get("value"))
radiusPlasmodesm_disp=float(etree.parse(dir + 'in/General.xml').getroot().xpath('radiusPlasmodesm_disp')[0].get("value"))

#Import Geometrical data
Plant=etree.parse(dir + 'in/' + Geom).getroot().xpath('Plant')[0].get("value")
path=etree.parse(dir + 'in/' + Geom).getroot().xpath('path')[0].get("value")
im_scale=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('im_scale')[0].get("value"))
Barrierrange=etree.parse(dir + 'in/' + Geom).getroot().xpath('Barrierrange/Barriers')
passage_cell_range=etree.parse(dir + 'in/' + Geom).getroot().xpath('passage_cell_range/passage_cell')
passage_cell_ID=[]
for passage_cell in passage_cell_range:
    passage_cell_ID.append(int(passage_cell.get("id")))
InterC_perim1=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('InterC_perim1')[0].get("value"))
InterC_perim2=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('InterC_perim2')[0].get("value"))
InterC_perim3=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('InterC_perim3')[0].get("value"))
InterC_perim4=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('InterC_perim4')[0].get("value"))
InterC_perim5=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('InterC_perim5')[0].get("value"))
kInterC=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('kInterC')[0].get("value"))
height=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('height')[0].get("value"))
thickness=float(etree.parse(dir + 'in/' + Geom).getroot().xpath('thickness')[0].get("value"))

#Import cellset data
tree = etree.parse(dir + 'cellsetdata/' + path) #Maize_Charles\\Maize_pip_cross4.xml') # #Parse literally decrypts the tree element data         SteleOK_high.xml
rootelt = tree.getroot()
waalls = rootelt.xpath('cells/cell/walls') #waalls contains cell wall groups info (one group by cell), searched by xpath ("Smart" element identifier)

#Set path
points = rootelt.xpath('walls/wall/points') #points contains the wall elements attributes
woual = rootelt.xpath('cells/cell/walls/wall') #woual contains the individual cell to wall associations
coul=rootelt.xpath('cells/cell') #coul contains the cell attributes
newpath=dir+Output_path
if not os.path.exists(newpath):
    os.makedirs(newpath)

#Coordinates
listee=[] #listeste=[]
for i in range(0, len(points)): #Loop on wall elements (groups of points)
	firstp = points[i].getchildren()[0] #Points to the first point of the wall element
	firstid = points[i].getparent().get("id") #Contains the wall id number
	#firstgp = points[i].getparent().get("group") #Contains the wall group number
	x1 = im_scale*int(firstp.get("x")) #X coordinate of the first point of the wall (micrometers)
	y1 = im_scale*int(firstp.get("y")) #Y coordinate of the first point of the wall (micrometers)
	#print("x1, y1, firstid", x1, y1, firstid)
	
	npo = len(points[i].getchildren()) #Total number of points in the ith wall element	
	lastp = points[i].getchildren()[npo - 1] #Points to the last point of the wall
	lastid = points[i].getparent().get("id") #Contains the wall id number
	#lastgp = points[i].getparent().get("group") #Contains the wall group number
	x2 = im_scale*int(lastp.get("x")) #X coordinate of the last point of the wall
	y2 = im_scale*int(lastp.get("y")) #Y coordinate of the last point of the wall
	#print("x2, y2, lastp", x2, y2, lastid)
	z1 = "x"+str(x1)+"y"+str(y1)
	z2 = "x"+str(x2)+"y"+str(y2)
	listee.append([z1, x1, y1, firstid])
	listee.append([z2, x2, y2, lastid])
	#listeste.append([z1, x1, y1, firstid, firstgp])
	#listeste.append([z2, x2, y2, lastid, lastgp])

G = nx.Graph() #Initializes graphical structure

#Add nodes
borlist1=[]
az=[]
xapp=[]
yapp=[]
length_tot=0.0 # (micrometers)
min_x_wall=inf
max_x_wall=0
for p in points: #Loop on wall elements (groups of points)
        pid= int((p.getparent().get)("id")) #pid records the current wall id number
        #print("wall id", pid, "number of points", len(p)) #,(p.getparent().get)("group")
        er=len(p) #Number of points in the wall element
        xsum=0.0
        ysum=0.0
        xprev=inf
        yprev=inf
        length=0.0
        for r in p: #Loop on points within the wall element to calculate their average X and Y coordinates 
                #print(pid, r.get("x"), r.get("y"))	#print("testX",(p.getchildren())("x"))
		x= im_scale*int(r.get("x")) #X coordinate of the point
		y= im_scale*int(r.get("y")) #X coordinate of the point
		xsum+= x
		ysum+= y
		if not xprev==inf:
		    length+=sqrt(square(x-xprev)+square(y-yprev))
		xprev=x
		yprev=y
        bx=xsum/er #Finalizing the averaging
        by=ysum/er
        length_tot+=length
        min_x_wall=min(min_x_wall,bx)
        max_x_wall=max(max_x_wall,bx)
	
        tuy=(p.getparent().get)("id"), bx, by #Wall id number and X Y position
        az.append(tuy) #az records walls id numbers and positions
        xapp.append(bx)#xapp records walls X positions
        yapp.append(by)#yapp records walls Y positions
	
        cooot=0
        for w in woual: #Loop on walls, by cell - wall association, hence a wall can be repeated if associated to two cells
            wid= int(w.get("id")) #Wall id number
            if wid == pid:
		cooot += 1 #Count cell - wall associations, for each wall
		#print("if wid = pid : coot+1", wid, pid,)
        if cooot == 1: #If we have only one cell associated to a wall element (soil-root interface or interface with xylem)
		borlin = 1
		#print("wid pid if coot = 1 : borkin 1", wid, pid,)
		borlist1.append(pid) #Lists the IDs of walls at soil-root interface or interface with xylem
        else: #If we have several cells associated to one wall element
		borlin = 0
		#print("wid pid else : borkin 0 et cooot", wid, pid, cooot) #Here wid is the last wall id due to the loop
        #print("wid, pid, cooot, borlin",wid, pid, cooot, borlin)
        #print("G vrai noeud",pid, borlin )		

        G.add_node(pid, indice=pid, borderlink=borlin, type="apo", position=(bx,by), height=height, length=length, thickness=thickness) #Saving wall attributes for graphical display (id, border, type, X and Y coordinates)

coro=int(G.number_of_nodes()) #Total number of walls included in G (not counting junction nodes)


dico={}
dico2={}
dico3={}
dico4={}

cr=coro
#Junction2Wall1=zeros((3*coro,1),dtype=int) #Provides one connected wall per junction
Junction2Wall1=[]
for a,b,c,d in listee: #a is the character string for X Y coordinates of the cell wall start or end, b is the X coordinate, c the Y coordinate, and d the cell wall ID number
        #print(a,b,c,d)
	if a in dico:
		dico[a].append(d) #Several cell wall ID numbers can correspond to the same X Y coordinate where they meet
	else:
		cr+=1 #New node count
		dico[a]=[d] #Saves the cell wall ID number associated to the X Y coordinate
		dico2[a]=[b,c] #Saves the X and Y coordinates of junctions the first time they appear in the loop
		bb=int(b) #X coordinate of the junction
		cc=int(c) #Y coordinate of the junction
		G.add_node(cr-1, indice=cr-1,borderlink=0, type="apo", position=(bb,cc), height=height, thickness=thickness) #Nodes are added at walls junctions (previous nodes corresponded to walls middle points). By default, borderlink is 0, but will be adjusted in next loop
		dico3[a]=[cr-1] #Saves the junction index
		#Junction2Wall1[cr-1]=d #We only need 1 connection to a wall by junction so that junctions at the interface with xylem could be connected to a wall that is connected to a cell that is associated to a cell layer
		Junction2Wall1.append([cr-1,int(d)])

ret=[]

for cle, valeur in dico.items(): #Loop on junctions
	dfg=0
	for x in valeur: #valeur is the list of cell walls ID meeting at the current junction
		if int(x) in borlist1: #If the cell wall has a junction at soil-root interface or interface with xylem
			dfg +=1 #Counting junctions at interface for current cell wall
	if dfg >= 2: #If the cell wall has both extremities at soil-root interface or at interface with xylem
		ret.append(cle) #Saving the X Y coordinates of junctions of walls delimiting interfaces


for cle, valeur in dico3.items(): #Loop on junctions
        #print(cle,valeur)
	if cle in ret: #cle contains the X Y coordinates of junctions, and ret the XY coordinates of junctions at interfaces
		for x in valeur: #valeur is the current junction ID, but has a list class, this loop take the single element of the list as single string that can be turned into an integer
		      G.node[int(x)]['borderlink'] = 1 #Interface junction nodes now have the value 1 in the attribute "borderlink"

zoro=int((G.number_of_nodes())) #Total number of nodes, including cell wall nodes and junction nodes


#Get X and Y for Cell nodes and cell nodes
for w in waalls: #Loop on cells. waalls contains cell wall groups info (one group by cell)
	totx=0.0 #Summing up cell walls X positions
	toty=0.0 #Summing up cell walls Y positions
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	cgroup=int(w.getparent().get("group")) #Cell type (1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
	div=float(len(w)) #Total number of walls around the current cell
	for r in w: #w points to the cell walls around the current cell
		#cellnumber2 = int(w.getparent().get("id")) #
		wid= int(r.get("id")) #Wall ID number
		#print (cellnumber1,cellnumber2, div, wid)
		#if cellnumber1== cellnumber2:
		totx += xapp[wid]#xapp contains the walls average X positions
		toty += yapp[wid]#yapp contains the walls average Y positions
	#print(totx, toty)
	finalx=totx/div #Average cell X position (from the average position of its walls)
	finaly=toty/div #Average cell Y position (from the average position of its walls)
	#print("cell:", cellnumber1, "x:", finalx, "y:", finaly)
	G.add_node(zoro + cellnumber1, indice=(zoro) + cellnumber1, borderlink=0, type="cell", position = (finalx,finaly), cgroup=cgroup) #Adding cell nodes
	#print((zoro) + int((w.getparent().get)("id")), cellnumber1)

Ncells=len(G.node)-zoro #Total number of cell walls
position=nx.get_node_attributes(G,'position') #Nodes XY positions (micrometers)

#add Edges
lat_dists=zeros((coro,1))
Nmb=0 #Total number of membranes
lengths=nx.get_node_attributes(G,'length') #Walls lengths (micrometers)
cellperimeter=np.linspace(0,0,Ncells)
for w in waalls: #Loop on cells. waalls contains cell wall groups info (one group by cell)
           	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
           	for r in w: #Loop for wall elements around the cell
          		wid= int(r.get("id")) #Cell wall ID
          		dist_cell=sqrt(square(position[wid][0]-position[zoro+cellnumber1][0])+square(position[wid][1]-position[zoro+cellnumber1][1])) #distance between wall node and cell node (micrometers)
          		lat_dists[wid]+=dist_cell
          		Nmb+=1
          		G.add_edge(zoro + cellnumber1, wid, path='membrane', height=height, length=lengths[wid], dist=dist_cell) #Adding all cell to wall connections (edges) #kaqp=kaqp, kw=kw, kmb=kmb, 
          		cellperimeter[cellnumber1]+=lengths[wid]

Cell_connec=-ones((Ncells,9),dtype=int) #Connected cells for further ranking
nCell_connec=zeros((Ncells,1),dtype=int) #Quantity of cell to cell connections
Wall2Cell1=zeros((zoro,1),dtype=int) #Provides one connected cell per wall node
for i in range(0, len(woual)): #Loop on walls, by cell - wall association, hence a wall can be repeated if associated to two cells. Parent structure: Cell/Walls/Wall
           	r1 = woual[i] #Points to the current wall
           	cellid1 = r1.getparent().getparent().get("id") #Cell1 ID number
           	id1 = r1.get("id") #Wall1 ID number
           	Wall2Cell1[int(id1)] = cellid1 #Walls at boundaries will send to the only cell they are associated with
           	for j in range(i + 1, len(woual) ): #Loop on cell-wall associations that are further down in the list
          		r2 = woual[j] #Points to the further down wall in the list of cell-wall associations
          		cellid2 = r2.getparent().getparent().get("id") #Cell2 ID number
          		id2 = r2.get("id") #Wall2 ID number
          		if id1 == id2: #If walls 1 and 2 are the same, then cells 1 and 2 are connected by plasmodesmata
         			#print(cellid1, cellid2, id1)
         			#print(zoro + int(cellid1),   zoro + int(cellid2))
         			G.add_edge(zoro + int(cellid1), zoro + int(cellid2), path='plasmodesmata', height=height, length=lengths[int(id1)]) #Adding all cell to cell connections (edges) #kpl=kpl, 
         			Cell_connec[int(cellid1)][nCell_connec[int(cellid1)]]=int(cellid2)
         			nCell_connec[int(cellid1)]+=1
         			Cell_connec[int(cellid2)][nCell_connec[int(cellid2)]]=int(cellid1)
         			nCell_connec[int(cellid2)]+=1

ret=[]

for junctionid, wallid in Junction2Wall1: #Associates a cell id to each junction so that we could associate junctions to cell layers eventually
    Wall2Cell1[int(junctionid)] = Wall2Cell1[int(wallid)]

for cle, valeur in dico.items(): #Loop on junctions between walls
        for cle3, valeur3 in dico3.items():
            if cle3 == cle: #Both cle and cle3 contain the X Y coordinates of wall junctions
                for x in valeur: #valeur is the list of cell walls ID meeting at the junction cle
                    G.add_edge(valeur3[0], int(x), path='wall', height=height, length=lengths[int(x)]/2, lat_dist=lat_dists[int(x)][0], thickness=thickness) #Adding junction to wall connections (edges) #kw=kw, 


#And calculation of the centre of gravity of the endodermis
x_grav=0.0 # (micrometers)
y_grav=0.0 # (micrometers)
n_cell_endo=0 #Counting the total number of cells in the endodermis
for w in waalls: #Loop on cells. waalls contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	if G.node[zoro + cellnumber1]['cgroup']==3: #Endodermis
	        x_grav+=position[zoro + cellnumber1][0]
	        y_grav+=position[zoro + cellnumber1][1]
	        n_cell_endo+=1
x_grav/=n_cell_endo
y_grav/=n_cell_endo

#Calculation of a and b parameter for cortex AQP activity radial distribution
Cell_rank=zeros((Ncells,1)) #Ranking of cells (1=Exodermis, 2=Epidermis, 3=Endodermis, 4*=Cortex, 5*=Stele, 16=Pericycle), stars are replaced by the ranking within cortical cells and stele cells
Layer_dist=zeros((62,1)) #Average cell layers distances from center of gravity, by cells ranking (1=Exodermis, 2=Epidermis, 3=Endodermis, 4*=Cortex, 5*=Stele, 16=Pericycle), stars are replaced by the ranking within cortical cells and stele cells
nLayer=zeros((62,1)) #Total number of cells in each rank (indices follow ranking numbers)
#angle_dist_endo_grav=array([-4,0]) #array of distances and angles between endo cells and grav. Initializing the array with values that will eventualy be deleted
#angle_dist_exo_grav=array([-4,0]) #array of distances and angles between exo cells and grav. Initializing the array with values that will eventualy be deleted
for w in waalls: #Loop on cells. waalls contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	celltype=G.node[zoro + cellnumber1]['cgroup'] #Cell type
	Cell_rank[cellnumber1]=celltype #Later on, cell types 4 and 5 will be updated to account for their ranking within the cortex / stele
	if celltype==16: #Pericycle
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[celltype]+=dist
	        nLayer[celltype]+=1
	elif celltype==3: #Endodermis
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        #angle=math.atan2(x_cell-x_grav,y_cell-y_grav) #radians between -pi and +pi
	        #angle_dist_endo_grav=vstack((angle_dist_endo_grav,array([angle,dist])))
	        Layer_dist[celltype]+=dist
	        nLayer[celltype]+=1
	elif celltype==1: #Exodermis
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        #angle=math.atan2(x_cell-x_grav,y_cell-y_grav) #radians between -pi and +pi
	        #angle_dist_exo_grav=vstack((angle_dist_exo_grav,array([angle,dist])))
	        Layer_dist[celltype]+=dist
	        nLayer[celltype]+=1
	elif celltype==2: #Epidermis
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[celltype]+=dist
	        nLayer[celltype]+=1

if nLayer[16]==0: #If there is no labelled pericycle
    stele_connec_rank=3 #Endodermis connected to stele cells
else:
    stele_connec_rank=16 #Pericycle connected to stele cells
if nLayer[1]==0: #If there is no labelled exodermis
    outercortex_connec_rank=2 #Cortex connected to epidermis cells
else:
    outercortex_connec_rank=1 #Cortex connected to exodermis cells
rank_cellperimeters_in=linspace(nan,nan,100)
rank_cellperimeters_out=linspace(nan,nan,100)
InterCid=list()
for w in waalls: #Loop on cells. waalls contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	celltype=Cell_rank[cellnumber1] #Cell types 4 and 5 updated to account for their ranking within the cortex / stele
	if celltype==4: #Cortex
	    if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==3): #Cell to cell connection with endodermis
	        Cell_rank[cellnumber1]=40
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[40]+=dist
	        nLayer[40]+=1
	        rank_cellperimeters_in[int(nLayer[40]-1)]=cellperimeter[cellnumber1]
	        if cellperimeter[cellnumber1]<InterC_perim1:
	            InterCid.append(cellnumber1)
	    elif any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==outercortex_connec_rank): #Cell to cell connection with exodermis
	        Cell_rank[cellnumber1]=49
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[49]+=dist
	        nLayer[49]+=1
	        rank_cellperimeters_out[int(nLayer[49]-1)]=cellperimeter[cellnumber1]
	        if cellperimeter[cellnumber1]<InterC_perim5:
	            InterCid.append(cellnumber1)
	elif celltype==5: #Stele
	    if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==stele_connec_rank): #Cell to cell connection with pericycle
	        Cell_rank[cellnumber1]=50
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[50]+=dist
	        nLayer[50]+=1

cortex_cellperimeters_in=rank_cellperimeters_in #Inner part of the cortex (close to endodermis)
cortex_cellperimeters_out=rank_cellperimeters_out #Outer part of cortex
for i in range(12):
    rank_cellperimeters_in=linspace(nan,nan,100)
    rank_cellperimeters_out=linspace(nan,nan,100)
    for w in waalls: #Loop on cells. waalls contains cell wall groups info (one group by cell)
        cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	celltype=Cell_rank[cellnumber1] #Cell types 4 and 5 updated to account for their ranking within the cortex / stele
	if celltype==4 and i<4: #Cortex
	   # if i<4: #Within 3 layers of cortical sides
	        if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==(40+i)): #Cell to cell connection with endodermis
	            Cell_rank[cellnumber1]=41+i
	            x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	            y_cell=position[zoro + cellnumber1][1]
	            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	            Layer_dist[41+i]+=dist
	            nLayer[41+i]+=1
	            rank_cellperimeters_in[int(nLayer[41+i]-1)]=cellperimeter[cellnumber1]
	            if i==0 and cellperimeter[cellnumber1]<InterC_perim2:
	                InterCid.append(cellnumber1)
	            elif i==1 and cellperimeter[cellnumber1]<InterC_perim3:
	                InterCid.append(cellnumber1)
	            elif i==2 and cellperimeter[cellnumber1]<InterC_perim4:
	                InterCid.append(cellnumber1)
	            elif i>2 and cellperimeter[cellnumber1]<InterC_perim5:
	                InterCid.append(cellnumber1)
	        elif any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==(49-i)): #Cell to cell connection with exodermis
	            Cell_rank[cellnumber1]=48-i
	            x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	            y_cell=position[zoro + cellnumber1][1]
	            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	            Layer_dist[48-i]+=dist
	            nLayer[48-i]+=1
	            rank_cellperimeters_out[int(nLayer[48-i]-1)]=cellperimeter[cellnumber1]
	            if cellperimeter[cellnumber1]<InterC_perim5:
	                InterCid.append(cellnumber1)
	   # else: #Center of cortex
	   #     Cell_rank[cellnumber1]=44
	   #     x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	   #     y_cell=position[zoro + cellnumber1][1]
	   #     dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	   #     Layer_dist[44]+=dist
	   #     nLayer[44]+=1
	elif celltype==5: #Stele
	    if i<10:
	        if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==(50+i)): #Cell to cell connection with pericycle
	            Cell_rank[cellnumber1]=51+i
	            x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	            y_cell=position[zoro + cellnumber1][1]
	            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	            Layer_dist[51+i]+=dist
	            nLayer[51+i]+=1
	    else: #No more than 11 stele cell layers
	        Cell_rank[cellnumber1]=61
	        x_cell=position[zoro + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[zoro + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[61]+=dist
	        nLayer[61]+=1
    if i<4:
        cortex_cellperimeters_in=vstack((cortex_cellperimeters_in,rank_cellperimeters_in))
        cortex_cellperimeters_out=vstack((rank_cellperimeters_out,cortex_cellperimeters_out))
cortex_cellperimeters=vstack((cortex_cellperimeters_in,cortex_cellperimeters_out))

#Calculating cell surfaces at tissue interfaces (total and interfacing with a cell that is not an intercellular space)
indice=nx.get_node_attributes(G,'indice') #Node indices (walls, junctions and cells)
Length_outer_cortex_tot=0.0 #Total cross-section membrane length at the interface between exodermis and cortex
Length_cortex_cortex_tot=0.0 #Total cross-section membrane length at the interface between cortex and cortex
Length_cortex_endo_tot=0.0 #Total cross-section membrane length at the interface between cortex and endodermis
Length_outer_cortex_nospace=0.0 #Cross-section membrane length at the interface between exodermis and cortex not including interfaces with intercellular spaces
Length_cortex_cortex_nospace=0.0 #Cross-section membrane length at the interface between exodermis and cortex not including interfaces with intercellular spaces
Length_cortex_endo_nospace=0.0 #Cross-section membrane length at the interface between exodermis and cortex not including interfaces with intercellular spaces
for node, edges in G.adjacency_iter() :
    i=indice[node] #Node ID number
    if i>=zoro: #Cell
        if G.node[i]['cgroup']==outercortex_connec_rank or G.node[i]['cgroup']==4 or G.node[i]['cgroup']==3: #exodermis or cortex or endodermis (or epidermis if there is no exodermis)
            if i-zoro not in InterCid: #The loop focuses on exo, cortex and endodermis cells that are not intercellular spaces
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == "plasmodesmata": #Plasmodesmata connection  #eattr is the edge attribute (i.e. connection type)
                        j = (indice[neighboor]) #neighbouring node number
                        l_membrane=eattr['length']
                        if (G.node[i]['cgroup']==outercortex_connec_rank and G.node[j]['cgroup']==4) or (G.node[j]['cgroup']==outercortex_connec_rank and G.node[i]['cgroup']==4):#Exodermis to cortex cell or vice versa (epidermis if no exodermis exists)
                            Length_outer_cortex_tot+=l_membrane
                            if j-zoro not in InterCid:
                                Length_outer_cortex_nospace+=l_membrane
                        elif (G.node[i]['cgroup']==4 and G.node[j]['cgroup']==4):#Cortex to cortex cell
                            Length_cortex_cortex_tot+=l_membrane
                            if j-zoro not in InterCid:
                                Length_cortex_cortex_nospace+=l_membrane
                        elif (G.node[i]['cgroup']==3 and G.node[j]['cgroup']==4) or (G.node[j]['cgroup']==3 and G.node[i]['cgroup']==4):#Cortex to endodermis cell or vice versa
                            Length_cortex_endo_tot+=l_membrane
                            if j-zoro not in InterCid:
                                Length_cortex_endo_nospace+=l_membrane

for i in range(62): #Finalizing distance averaging
    if nLayer[i]>0:
        Layer_dist[i]=Layer_dist[i]/nLayer[i]

#Discretization based on effective cell layering
r_discret=array([0])
j=0 #Counts cell layers
k=0 #Counts tissue types
rank2row=zeros((62,1))*nan
Layer_dist2=array([0])
for i in range(61, 49, -1): #rank2row list, Stele
    if nLayer[i]>0:
        rank2row[i]=j
        Layer_dist2=vstack((Layer_dist2,Layer_dist[i]))
        j+=1
k+=1
r_discret=vstack((r_discret,j))
if nLayer[16]>0:
    rank2row[16]=j #Pericycle
    Layer_dist2=vstack((Layer_dist2,Layer_dist[16]))
    j+=1
    k+=1
    r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
rank2row[3]=j #2 rows for endodermis (inner and outer) + 2 rows for passage cells in between
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
j+=4
k+=1
r_discret=vstack((r_discret,j-sum(r_discret[1:k])))

i1=40
nLayer_ref=nLayer[i1]
Layer_dist_ref=Layer_dist[i1]
rank2row[i1]=j
Layer_dist2=vstack((Layer_dist2,Layer_dist[i1]))
j+=1 #number of the row that was just added
i1=41
ratio_complete=0.75
while i1<50: #Cortex
    if nLayer[i1]>ratio_complete*nLayer_ref: #Likely complete rank/layer
        rank2row[i1]=j
        Layer_dist2=vstack((Layer_dist2,Layer_dist[i1]))
        j+=1
        nLayer_ref=nLayer[i1]
        Layer_dist_ref=Layer_dist[i1]
        i1+=1
    elif nLayer[i1]>0: #Likely incomplete rank/layer i1
        #Find next non-empty rank (i1+i2)
        i2=i1+1
        while nLayer[i2]==0:
            i2+=1
        #Check if incomplete layers/ranks i1 and i2 would constitute a full layer
        if nLayer[i1]+nLayer[i2]>ratio_complete*nLayer_ref: #They do
            #Does layer/rank i2 constitute a full layer?
            if nLayer[i2]>ratio_complete*nLayer_ref: #Yes
                #Rank i1 is added to the closest layer
                if abs(Layer_dist[i1]-Layer_dist[i2])<abs(Layer_dist[i1]-Layer_dist_ref): #i2 is the closest layer
                    rank2row[i1]=j
                    rank2row[i2]=j
                    Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2])/(nLayer[i1]+nLayer[i2])
                    Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                    j+=1
                else: #i1-1 is the closest layer to rank i1
                    rank2row[i1]=j-1
                    Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer_ref*Layer_dist_ref)/(nLayer[i1]+nLayer_ref)
                    Layer_dist2[j-1]=Layer_dist_avg
                    rank2row[i2]=j
                    Layer_dist2=vstack((Layer_dist2,Layer_dist[i2]))
                    j+=1
            else: #Rank i2 does not constitute a full layer, then i1 and i2 form a layer together
                rank2row[i1]=j
                rank2row[i2]=j
                Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2])/(nLayer[i1]+nLayer[i2])
                Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                j+=1
            i1=i2+1
        else: #Ranks i1 and i2 don't constitute a full layer, then i1 and i2 can either be both added to the same layer, or to separate layers
            #Ranks i1 and i2 are added to the closest layer
            if nLayer[i2+1]>0 and i2+1<50:#the next full layer needs to be part of the cortex
                if abs(Layer_dist[i1]-Layer_dist[i2+1])<abs(Layer_dist[i1]-Layer_dist_ref): #i2+1 is the closest layer, so i1 and i2 are both added with i2+1
                    rank2row[i1]=j
                    rank2row[i2]=j
                    rank2row[i2+1]=j
                    Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2]+nLayer[i2+1]*Layer_dist[i2+1])/(nLayer[i1]+nLayer[i2]+nLayer[i2+1])
                    Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                    j+=1
                else:#i1 is closer to i1-1 than i2+1
                    rank2row[i1]=j #Correction 18/04: j replaces j-1
                    if abs(Layer_dist[i2]-Layer_dist[i2+1])<abs(Layer_dist[i2]-Layer_dist_ref): #i2 is closer to i2+1 than to i1-1
                        Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer_ref*Layer_dist_ref)/(nLayer[i1]+nLayer_ref)
                        Layer_dist2[j]=Layer_dist_avg #Correction 18/04: j replaces j-1
                        rank2row[i2]=j
                        rank2row[i2+1]=j
                        Layer_dist_avg=(nLayer[i2]*Layer_dist[i2]+nLayer[i2+1]*Layer_dist[i2+1])/(nLayer[i2]+nLayer[i2+1])
                        Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                        j+=1
                    else: #i2 is closer to i1-1 than to i2+1
                        rank2row[i2]=j-1
                        Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2]+nLayer_ref*Layer_dist_ref)/(nLayer[i1]+nLayer[i2]+nLayer_ref)
                        Layer_dist2[j]=Layer_dist_avg #Correction 18/04: j replaces j-1
                        rank2row[i2+1]=j
                        Layer_dist2=vstack((Layer_dist2,Layer_dist[i2+1]))
                        j+=1
            else: #We just merge i1 and i2
                rank2row[i1]=j
                rank2row[i2]=j
                Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2])/(nLayer[i1]+nLayer[i2])
                Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                j+=1
            i1=i2+2
    elif nLayer[i1]==0: #No more cortex layer (likely there was only one of them)
        i1=50
#for i in range(40, 50, 1): #Cortex
#    if nLayer[i]>0:
#        rank2row[i]=j
#        Layer_dist2=vstack((Layer_dist2,Layer_dist[i]))
#        j+=1
k+=1
r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
if nLayer[1]>0: #if there is an exodermal layer
    rank2row[1]=j #2 rows for exodermis (inner and outer)
    Layer_dist2=vstack((Layer_dist2,Layer_dist[1]))
    Layer_dist2=vstack((Layer_dist2,Layer_dist[1]))
    j+=2
    k+=1
    r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
rank2row[2]=j #Epidermis
Layer_dist2=vstack((Layer_dist2,Layer_dist[2]))
j+=1
k+=1
r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
r_discret[0]=j
Layer_dist2=Layer_dist2[1:]


#Then min and max distances to cross-section centre
dmax_cortex=0.0 #For gradient of relative AQP distribution
dmin_cortex=inf #For gradient of relative AQP distribution
davg_epi=0.0 #avg distance of between exodermis and centre of gravity (micrometers)
dist_grav=zeros((coro,1)) #distance between node and cross-section center of gravity (micrometers)
temp=0.0 #Counting the number of membranes from epidermis
for w in waalls: #Loop on cells. waalls contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	for r in w: #Loop for wall elements around the cell
		wid= int(r.get("id")) #Cell wall ID
		#angle=math.atan2(position[wid][0]-x_grav,position[wid][1]-y_grav) #radians between -pi and +pi
		dist_grav[wid]=sqrt(square(position[wid][0]-x_grav)+square(position[wid][1]-y_grav)) #Distance between membrane and cross section centre of gravity (micrometers)
		if G.node[zoro + cellnumber1]['cgroup']==4: #Cortex
		    dmax_cortex=max(dmax_cortex,dist_grav[wid])
		    dmin_cortex=min(dmin_cortex,dist_grav[wid])
		elif G.node[zoro + cellnumber1]['cgroup']==2: #Epidermis
		    davg_epi+=dist_grav[wid]
		    temp+=1.0
davg_epi/=temp #Last step of averaging (note that we take both inner and outer membranes into account in the averaging)
perimeter=2*pi*davg_epi*1.0E-04 #(cm)

#Creating link to borders (soil and xylem)
borderlink=nx.get_node_attributes(G,'borderlink') #Is 1 at interfaces with soil and xylem for wall and junction nodes

#Calculating cell borders accounting for wall thickness
ThickWalls=[]
nThickWalls=zeros((2*coro,1)) #This will save how many new junction wall ID were already saved in connections to new wall ID (the vector is a little too long)
ThickWallsX=[]
Wall2NewWall=empty((coro,2))
Wall2NewWall[:]=NAN
nWall2NewWall=zeros((coro,1))
Wall2NewWallX=empty((cr,8)) #This one also includes "junction" to "new junctions" ID
Wall2NewWallX[:]=NAN
nWall2NewWallX=zeros((cr,1))
ThickWallPolygonX=empty((coro*2,4))
ThickWallPolygonX[:]=NAN
nThickWallPolygonX=zeros((coro*2,1))
Wall2Cell=empty((coro,2))
Wall2Cell[:]=NAN
nWall2Cell=zeros((coro,1))
Junction2Wall2Cell=empty((cr,12))
Junction2Wall2Cell[:]=NAN
nJunction2Wall2Cell=zeros((cr,1))
Junction2Wall=empty((cr,12))
Junction2Wall[:]=NAN
nJunction2Wall=zeros((cr,1))
Wall2Junction=empty((coro,2))
Wall2Junction[:]=NAN
nWall2Junction=zeros((coro,1))
Cell2ThickWalls=empty((Ncells,20))
Cell2ThickWalls[:]=NAN
nCell2ThickWalls=zeros((Ncells,1))
twpid=0 #Thick wall point ID
twpidX=0
for node, edges in G.adjacency_iter():
        i=indice[node]
        if i<coro: #wall that is not a junction (connected to a cell)
            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                j=indice[neighboor]
                if G.node[j]['type']=='cell':
                    if not position[j][0]==position[i][0]:
                        slopeCG=(position[j][1]-position[i][1])/(position[j][0]-position[i][0]) #slope of the line connecting the wall node to the center of gravity of the cell
                    else:
                        slopeCG=inf
                    x=position[i][0]+cos(arctan(slopeCG))*thickness_disp/2*sign(position[j][0]-position[i][0]) #-> new position of the wall node on the current cell side
                    y=position[i][1]+sin(arctan(slopeCG))*thickness_disp/2*sign(position[j][0]-position[i][0])
                    ThickWalls.append(array((twpid,i,j,x,y,inf,inf,G.node[i]['borderlink']))) #Adds the thick wall node ID, its parent node ID, the associated cell ID, the new coordinates in X and Y, and the two neighbouring new junction walls IDs (not known at this point in the loop)
                    Cell2ThickWalls[j-zoro][int(nCell2ThickWalls[j-zoro])]=twpid
                    nCell2ThickWalls[j-zoro]+=1
                    ThickWallsX.append((twpidX,x,y,i,j)) #Adds the thick wall node ID, the new coordinates in X and Y, its parent node ID, the associated cell ID
                    ThickWallPolygonX[2*i][int(nThickWallPolygonX[2*i])]=twpidX
                    ThickWallPolygonX[2*i+1][int(nThickWallPolygonX[2*i+1])]=twpidX
                    nThickWallPolygonX[2*i]+=1
                    nThickWallPolygonX[2*i+1]+=1
                    Wall2Cell[i][int(nWall2Cell[i])]=j
                    nWall2Cell[i]+=1
                    Wall2NewWall[i][int(nWall2NewWall[i])]=twpid
                    nWall2NewWall[i]+=1  #the count in nWall2NewWall is actually the same as nWall2Cell 
                    Wall2NewWallX[i][int(nWall2NewWallX[i])]=twpidX
                    nWall2NewWallX[i]+=1
                    twpid+=1
                    twpidX+=1
                    if G.node[i]['borderlink']==1:
                        x=position[i][0]-cos(arctan(slopeCG))*thickness_disp/2*sign(position[j][0]-position[i][0]) #-> new position of the wall node opposite to the current cell side
                        y=position[i][1]-sin(arctan(slopeCG))*thickness_disp/2*sign(position[j][0]-position[i][0])
                        ThickWallsX.append((twpidX,x,y,i,inf))
                        ThickWallPolygonX[2*i][int(nThickWallPolygonX[2*i])]=twpidX
                        ThickWallPolygonX[2*i+1][int(nThickWallPolygonX[2*i+1])]=twpidX
                        nThickWallPolygonX[2*i]+=1
                        nThickWallPolygonX[2*i+1]+=1
                        Wall2NewWallX[i][int(nWall2NewWallX[i])]=twpidX
                        nWall2NewWallX[i]+=1
                        twpidX+=1
                elif G.node[j]['type']=='apo': #Node j is a junction
                    Wall2Junction[i][int(nWall2Junction[i])]=j
                    nWall2Junction[i]+=1
for node, edges in G.adjacency_iter():
        i=indice[node]
        if i<coro: #wall that is not a junction (connected to a cell)
            for neighboor, eattr in edges.items(): #Loop on connections (edges) separated to make sure that Wall2Cell[i] is complete
                j=indice[neighboor]
                if G.node[j]['type']=='apo': #then j is a junction node 
                    for cid in Wall2Cell[i]: #Wall2Cell[i][k] are the cell node ID associated to Wall i
                        if cid not in Junction2Wall2Cell[j]:
                            Junction2Wall2Cell[j][int(nJunction2Wall2Cell[j])]=cid #Writes the cells indirectly associated to each junction
                            nJunction2Wall2Cell[j]+=1
                            Junction2Wall[j][int(nJunction2Wall[j])]=i #Writes the walls directly associated to each junction
                            nJunction2Wall[j]+=1  #the count in nJunction2Wall is actually the same as nJunction2Wall2Cell 
                        else: #Cell already associated to junction j through another wall => we can consider that the junction j is directly associated to the cell Wall2Cell[i][k] geometrically
                            for id1, val in enumerate(Junction2Wall2Cell[j]): 
                                if val==cid: #Finding the position of the current cell in the list of cells related to the junction 
                                    wid1=int(Junction2Wall[j][id1]) #At the same position in Junction2Wall we can find the "other" wall ID that was already associated to the same junction and cell
                                    break
                            for id1, val in enumerate(Wall2Cell[wid1]): 
                                if val==cid: #Finding the position of the current cell in the list of cells related to the "other wall"
                                    twpid1=int(Wall2NewWall[wid1][id1]) #At the same position in Wall2NewWall we can find the "new wall" ID that was already associated to the same "other wall" and cell
                                    break
                            for id1, val in enumerate(Wall2Cell[i]):
                                if val==cid: #Finding the position of the current cell in the list of cells related to the "current wall"
                                    twpid2=int(Wall2NewWall[i][id1]) #At the same position in Wall2NewWall we can find the "new wall" ID that was already associated to the same "current wall" and cell
                                    break
                            if not position[cid][0]==position[j][0]:
                                slopeCG=(position[cid][1]-position[j][1])/(position[cid][0]-position[j][0]) #slope of the line connecting the junction node j to the center of gravity of the cell Wall2Cell[i][k]
                            else:
                                slopeCG=inf
                            x=position[j][0]+cos(arctan(slopeCG))*thicknessJunction_disp/2*sign(position[cid][0]-position[j][0])
                            y=position[j][1]+sin(arctan(slopeCG))*thicknessJunction_disp/2*sign(position[cid][0]-position[j][0])
                            ThickWalls.append(array((twpid,j,int(cid),x,y,twpid1,twpid2,G.node[j]['borderlink']))) #Adds the thick wall node ID, its parent node ID, the associated cell ID, the new coordinates in X and Y, and in this case, the two neighbouring walls that will consitute 2 neighbouring "cells" in the sense of pvtk
                            ThickWalls[twpid1][int(5+nThickWalls[twpid1])]=twpid
                            ThickWalls[twpid2][int(5+nThickWalls[twpid2])]=twpid
                            nThickWalls[twpid1]+=1
                            nThickWalls[twpid2]+=1
                            Cell2ThickWalls[int(cid-zoro)][int(nCell2ThickWalls[int(cid-zoro)])]=twpid
                            nCell2ThickWalls[int(cid-zoro)]+=1
                            ThickWallsX.append((twpidX,x,y,j,cid,i,wid1)) #Adds the thick wall node ID, the new coordinates in X and Y, its parent (junction) node ID, the associated cell ID, the two original neighbouring walls
                            #Which of the 2 polygons associated to wall node "i" do we add the twpid to?
                            #Only one of the is related to junction node "j". Need to look at position into Wall2Junction
                            for id1, val in enumerate(Wall2Junction[i]): 
                                if val==j: #Finding the position "id1" (index equal to 0 or 1) of the current junction in the list of junctions related to the wall "i" 
                                    break #2*i+id1 is the row of the polygon connecting wall "i" to junction "j"
                            if ThickWallsX[int(ThickWallPolygonX[int(2*i+id1)][1])][4]==cid: #The polygon nodes need to be consecutive ("no zigzag"), hence the ordering of the 3rd and 4th twpid matter
                                ThickWallPolygonX[int(2*i+id1)][2]=twpidX
                            else:
                                ThickWallPolygonX[int(2*i+id1)][3]=twpidX
                            nThickWallPolygonX[int(2*i+id1)]+=1
                            #We do the same for the other polygon that includes twpid (connects wall "wid1" to node "j")
                            for id1, val in enumerate(Wall2Junction[wid1]): 
                                if val==j: #Finding the position "id1" (index equal to 0 or 1) of the current junction in the list of junctions related to the wall "i" 
                                    break #2*i+id1 is the row of the polygon connecting wall "i" to junction "j"
                            if ThickWallsX[int(ThickWallPolygonX[int(2*wid1+id1)][1])][4]==cid: #The polygon nodes need to be consecutive ("no zigzag"), hence the ordering of the 3rd and 4th twpid matter
                                ThickWallPolygonX[int(2*wid1+id1)][2]=twpidX
                            else:
                                ThickWallPolygonX[int(2*wid1+id1)][3]=twpidX
                            nThickWallPolygonX[int(2*wid1+id1)]+=1
                            #New walls ID for each original wall
                            Wall2NewWallX[j][int(nWall2NewWallX[j])]=twpidX
                            nWall2NewWallX[j]+=1
                            #Increase the wall count
                            twpid+=1
                            twpidX+=1
                            if G.node[j]['borderlink']==1: #If the wall is at a border, there is only one cell => need to add the opposite new wall node independently
                                if not slopeCG==0.0:
                                    slopeCG=-1/slopeCG
                                else:
                                    slopeCG=inf
                                x=position[j][0]-cos(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[j][0]) #-> new position of the wall node opposite to the current cell side
                                y=position[j][1]-sin(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[j][0])
                                ThickWallsX.append((twpidX,x,y,j,inf,i))
                                #This new wall node is added to the border polygon. Is the border polygon associated to wall node "i" or "wid1"?
                                if G.node[i]['borderlink']==1:
                                    wid=i
                                else:
                                    wid=wid1
                                #Which of the 2 polygons associated to wall node "wid" do we add the twpid to?
                                #Only one of the is related to junction node "j". Need to look at position into Wall2Junction
                                for id1, val in enumerate(Wall2Junction[wid]): 
                                    if val==j: #Finding the position "id1" (index equal to 0 or 1) of the current junction in the list of junctions related to the wall "i" 
                                        break #2*i+id1 is the row of the polygon connecting wall "i" to junction "j"
                                if ThickWallsX[int(ThickWallPolygonX[int(2*wid+id1)][1])][4]==cid: #The polygon nodes need to be consecutive ("no zigzag"), hence the ordering of the 3rd and 4th twpid matter
                                    ThickWallPolygonX[int(2*wid+id1)][3]=twpidX
                                else:
                                    ThickWallPolygonX[int(2*wid+id1)][2]=twpidX
                                nThickWallPolygonX[int(2*wid+id1)]+=1
                                #New walls ID for each original wall
                                Wall2NewWallX[j][int(nWall2NewWallX[j])]=twpidX
                                nWall2NewWallX[j]+=1
                                twpidX+=1








#######################################
##Variables for potential calculation##
#######################################

#Import Hydraulic data
kwrange=etree.parse(dir + 'in/' + Hydr).getroot().xpath('kwrange/kw')
kw_barrier_range=etree.parse(dir + 'in/' + Hydr).getroot().xpath('kw_barrier_range/kw_barrier')
kmb=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('km')[0].get("value"))
kAQPrange=etree.parse(dir + 'in/' + Hydr).getroot().xpath('kAQPrange/kAQP')
Kplrange=etree.parse(dir + 'in/' + Hydr).getroot().xpath('Kplrange/Kpl')
Fplxheight=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight')[0].get("value"))
Fplxheight_epi_exo=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_epi_exo')[0].get("value"))
Fplxheight_outer_cortex=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_outer_cortex')[0].get("value"))*Length_outer_cortex_tot /Length_outer_cortex_nospace
Fplxheight_cortex_cortex=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_cortex_cortex')[0].get("value"))*Length_cortex_cortex_tot /Length_cortex_cortex_nospace
Fplxheight_cortex_endo=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_cortex_endo')[0].get("value"))*Length_cortex_endo_tot /Length_cortex_endo_nospace
Fplxheight_endo_endo=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_endo_endo')[0].get("value"))
Fplxheight_endo_peri=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_endo_peri')[0].get("value"))
Fplxheight_peri_peri=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_peri_peri')[0].get("value"))
Fplxheight_peri_stele=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_peri_stele')[0].get("value"))
Fplxheight_stele_stele=float(etree.parse(dir + 'in/' + Hydr).getroot().xpath('Fplxheight_stele_stele')[0].get("value"))
path_hydraulics=etree.parse(dir + 'in/' + Hydr).getroot().xpath('path_hydraulics/Output')

#Import boundary conditions
Psi_soil_range=etree.parse(dir + 'in/' + BC).getroot().xpath('Psi_soil_range/Psi_soil')
Psi_xyl_range=etree.parse(dir + 'in/' + BC).getroot().xpath('Psi_xyl_range/Psi_xyl')
Psi_cell_range=etree.parse(dir + 'in/' + BC).getroot().xpath('Psi_cell_range/Psi_cell')

#Unit changes
sperd=24.0*3600.0 #(seconds per day)
cmperm=100.0 #(cm per metre)

#Start the loop of hydraulic properties
Nhydraulics=len(kwrange) #Including the "forced scenario"
for h in range(Nhydraulics):
    newpath=dir+Output_path+path_hydraulics[h].get("path")
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    #System solving
    Nscenarios=len(Psi_soil_range) #Including the "forced scenario"
    Psi_soil=zeros((1,Nscenarios))
    Os_soil=zeros((2,Nscenarios))
    Psi_xyl=zeros((1,Nscenarios))
    Os_xyl=zeros((1,Nscenarios))
    Os_cortex=zeros((1,Nscenarios))
    s_factor=zeros((1,Nscenarios))
    Nbarriers=len(Barrierrange)
    UptakeLayer_plus=zeros((int(r_discret[0]),Nbarriers,Nscenarios))
    UptakeLayer_minus=zeros((int(r_discret[0]),Nbarriers,Nscenarios))
    Q_xyl_layer=zeros((int(r_discret[0]),Nbarriers,Nscenarios))
    STFcell_plus=zeros((Ncells,Nbarriers))
    STFcell_minus=zeros((Ncells,Nbarriers))
    STFlayer_plus=zeros((int(r_discret[0]),Nbarriers))
    STFlayer_minus=zeros((int(r_discret[0]),Nbarriers))
    PsiCellLayer=zeros((int(r_discret[0]),Nbarriers,Nscenarios))
    PsiWallLayer=zeros((int(r_discret[0]),Nbarriers,Nscenarios))
    OsCellLayer=zeros((int(r_discret[0]),Nscenarios))
    OsWallLayer=zeros((int(r_discret[0]),Nscenarios))
    NWallLayer=zeros((int(r_discret[0]),Nbarriers,Nscenarios))
    #UptakeDistri_plus=zeros((40,3,8))#the size will be adjusted, but won't be more than 40. Dimension 1: radial position, 2: compartment, 3: scenario
    #UptakeDistri_minus=zeros((40,3,8))
    Q_tot=zeros((Nbarriers,Nscenarios)) #(cm^3/d)
    kr_tot=zeros((Nbarriers,1))
    
    iBarriers=-1 #Iteration index in the Barriers loop
    for Barriers in Barrierrange: #Apoplastic barriers (1:Endodermis radial walls, 2:Endodermis with passage cells, 3: Endodermis full, 4: Endodermis full and exodermis radial walls)
        Maturity=int(Barriers.get("value"))
        
        #Index for barriers loop
        iBarriers+=1
        
        #Scenarios concern boundary conditions only
        count=0
        
        #Soil and xylem water potentials
        Psi_soil[0][count]=float(Psi_soil_range[count].get("pressure")) #Soil water potential (hPa)
        Psi_xyl[0][count]=float(Psi_xyl_range[count].get("pressure")) #Xylem water potential (hPa)
                
        #Cell wall hydraulic conductivity
        kw = float(kwrange[h].get("value"))
        kw_barrier = float(kw_barrier_range[h].get("value"))
        #kw_barrier = kw/10.0
        if Maturity==0: #No Casparian strip ###Yet to come: Punctured Casparian strip as in Steudle et al. (1993)
            kw_endo_endo=kw
            kw_puncture=kw
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_endo_peri=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Maturity==1: #Endodermis radial walls
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_endo_peri=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Maturity==2: #Endodermis with passage cells
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_endo_peri=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Maturity==3: #Endodermis full
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_endo_peri=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Maturity==4: #Endodermis full and exodermis radial walls
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_endo_peri=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        
        #Plasmodesmatal hydraulic conductance
        Kpl = float(Kplrange[h].get("value"))
        
        #Contribution of aquaporins to membrane hydraulic conductivity
        kaqp_cortex = float(kAQPrange[h].get("value"))
        
        ######################
        ##Filling the matrix##
        ######################
        
        matrix = np.zeros(((G.number_of_nodes()),G.number_of_nodes())) #Initializes the Doussan matrix
        
        Kmb=zeros((Nmb,1)) #Stores membranes conductances for the second K loop
        jmb=0 #Index of membrane in Kmb
        #c = csv.writer(open(dir + "matrice.csv", "wb")) #Creates a csv file to write the matrix
        #c.writerow (["i","j"])
        #c.writerow([[i], [j]])
        #Adding matrix components at cell-cell, cell-wall, and wall-junction connections
        for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
            i=indice[node] #Node ID number
            count_endo=0 #total number of endodermis cells around the wall
            count_peri=0 #total number of pericycle cells around the wall
            count_exo=0 #total number of exodermis cells around the wall
            count_epi=0 #total number of epidermis cells around the wall
            count_cortex=0 #total number of cortical cells around the wall
            count_passage=0 #total number of passage cells around the wall
            if i<coro: #wall ID
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == 'membrane': #Wall connection
                        if any(passage_cell_ID==array((indice[neighboor])-zoro)):
                            count_passage+=1
                        if G.node[neighboor]['cgroup']==3:#Endodermis
                            count_endo+=1
                        elif G.node[neighboor]['cgroup']==16 or G.node[neighboor]['cgroup']==5:#Pericycle or stele
                            count_peri+=1
                        elif G.node[neighboor]['cgroup']==4:#Cortex
                            count_cortex+=1
                        elif G.node[neighboor]['cgroup']==1:#Exodermis
                            count_exo+=1
                        elif G.node[neighboor]['cgroup']==2:#Epidermis
                            count_epi+=1
            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                j = (indice[neighboor]) #neighbouring node number
                if j > i: #Only treating the information one way to save time
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if path == 'wall': #Wall connection
                        #K = eattr['kw']*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                        if count_endo>=2: #wall between two endodermis cells
                            K = kw_endo_endo*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)  #(eattr['height']*eattr['thickness'])/eattr['length']#
                        elif count_peri>0 and count_endo>0: #wall between endodermis and pericycle
                            if count_passage>0:
                                K = kw_passage*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #(eattr['height']*eattr['thickness'])/eattr['length']#
                            else:
                                K = kw_endo_peri*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d) #(eattr['height']*eattr['thickness'])/eattr['length']#
                        elif count_peri==0 and count_endo==1: #wall between endodermis and cortex
                            if count_passage>0:
                                K = kw_passage*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length']  #(eattr['height']*eattr['thickness'])/eattr['length']#
                            else:
                                K = kw_endo_cortex*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)  #(eattr['height']*eattr['thickness'])/eattr['length']#
                        elif count_exo>=2: #wall between two exodermis cells
                            K = kw_exo_exo*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)  #(eattr['height']*eattr['thickness'])/eattr['length']#
                        else: #other walls
                            K = kw*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)  #(eattr['height']*eattr['thickness'])/eattr['length']#
                    elif path == "membrane": #Membrane connection
                        #K = (eattr['kmb']+eattr['kaqp'])*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                        if j-zoro in InterCid: #the neighbour is an intercellular space "cell". Between j and i connected by a membrane, only j can be cell because j>i
                            kaqp=kInterC
                        else:
                            kaqp=kaqp_cortex
                        if count_endo>=2: #wall between two endodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                            if kw_endo_endo==0.00:
                                K=0.00
                            else:
                                K = 1/(1/(kw_endo_endo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                        elif count_exo>=2: #wall between two exodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                            if kw_exo_exo==0.00:
                                K=0.00
                            else:
                                K = 1/(1/(kw_exo_exo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                        elif count_peri>0 and count_endo>0: #wall between endodermis and pericycle, in this case the suberized wall can limit the transfer of water between cell and wall
                            if count_passage>0:
                                K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                            else:
                                if kw_endo_peri==0.00:
                                    K=0.00
                                else:
                                    K = 1/(1/(kw_endo_peri/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                        elif count_peri==0 and count_endo==1: #wall between cortex and endodermis, in this case the suberized wall can limit the transfer of water between cell and wall
                            if kaqp==0.0:
                                K=1.00E-14
                            else:
                                if count_passage>0:
                                    K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                                else:
                                    if kw_endo_cortex==0.00:
                                        K=0.00
                                    else:
                                        K = 1/(1/(kw_endo_cortex/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                        else:
                            if kaqp==0.0:
                                K=1.00E-14
                            else:
                                K = 1/(1/(kw/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                        Kmb[jmb]=K
                        jmb+=1
                    elif path == "plasmodesmata": #Plasmodesmata connection
                        if j-zoro in InterCid or i-zoro in InterCid: #one of the connected cells is an intercellular space "cell".
                            K=0.0
                        elif (G.node[i]['cgroup']==2 and G.node[j]['cgroup']==1) or (G.node[j]['cgroup']==2 and G.node[i]['cgroup']==1):#Epidermis to exodermis cell or vice versa
                            K = Kpl*Fplxheight_epi_exo*1.0E-04*eattr['length'] #Will not be used in case there is no exodermal layer
                        elif (G.node[i]['cgroup']==outercortex_connec_rank and G.node[j]['cgroup']==4) or (G.node[j]['cgroup']==outercortex_connec_rank and G.node[i]['cgroup']==4):#Exodermis to cortex cell or vice versa
                            K = Kpl*Fplxheight_outer_cortex*1.0E-04*eattr['length']
                        elif (G.node[i]['cgroup']==4 and G.node[j]['cgroup']==4):#Cortex to cortex cell
                            K = Kpl*Fplxheight_cortex_cortex*1.0E-04*eattr['length']
                        elif (G.node[i]['cgroup']==3 and G.node[j]['cgroup']==4) or (G.node[j]['cgroup']==3 and G.node[i]['cgroup']==4):#Cortex to endodermis cell or vice versa
                            K = Kpl*Fplxheight_cortex_endo*1.0E-04*eattr['length']
                        elif (G.node[i]['cgroup']==3 and G.node[j]['cgroup']==3):#Endodermis to endodermis cell
                            K = Kpl*Fplxheight_endo_endo*1.0E-04*eattr['length']
                        elif (G.node[i]['cgroup']==3 and G.node[j]['cgroup']==16) or (G.node[j]['cgroup']==3 and G.node[i]['cgroup']==16):#Pericycle to endodermis cell or vice versa
                            K = Kpl*Fplxheight_endo_peri*1.0E-04*eattr['length']
                        elif (G.node[i]['cgroup']==16 and G.node[j]['cgroup']==5) or (G.node[j]['cgroup']==16 and G.node[i]['cgroup']==5):#Pericycle to stele cell or vice versa
                            K = Kpl*Fplxheight_peri_stele*1.0E-04*eattr['length']
                        elif (G.node[i]['cgroup']==5 and G.node[j]['cgroup']==5):#Stele to stele cell
                            K = Kpl*Fplxheight_stele_stele*1.0E-04*eattr['length']
                        else: #Default plasmodesmatal frequency
                            K = Kpl*Fplxheight*1.0E-04*eattr['length'] #eattr['kpl']
                    matrix[i][i] -= K #Filling the Doussan matrix (symmetric)
                    matrix[i][j] += K
                    matrix[j][i] += K
                    matrix[j][j] -= K
    
    
        #Adding matrix components at soil-wall and wall-xylem connections & rhs terms
        rhs_o = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of osmotic potentials
        rhs_p = np.zeros((G.number_of_nodes(),2)) #Initializing the right-hand side matrix of hydrostatic potentials
        #rhs_p = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of hydrostatic potentials
        listesol=[]
        listxy=[]
        K_soil=[]
        K_xyl=[]
        temp=0
        for w in woual: #Loop on cell-wall associations (cells/cell/walls/wall) necessary to know the cell layer (epidermis/endodermis)
            pgroup=int(w.getparent().getparent().get("group")) #Cell type (2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
            idwall=int(w.get("id")) #Wall ID number
            if G.node[idwall]['borderlink'] == 1:
                if pgroup == 2: #Wall (not including junctions) connected to soil
                    K=kw*1.0E-04*(G.node[idwall]['length']/2*G.node[idwall]['height'])/G.node[idwall]['thickness']
                    matrix[idwall][idwall] -= K #Doussan matrix
                    rhs_p[idwall][0] = -K    #Right-hand side vector, could become Psi_soil[idwall], which could be a function of the horizontal position
                    #rhs_p[idwall] -= Psi_soil[0][count]*K    #Right-hand side vector, could become Psi_soil[idwall], which could be a function of the horizontal position
                    listesol.append(idwall)
                    K_soil.append(K)
                elif pgroup == 5: #Wall (not including junctions) connected to xylem
                    K=kw*1.0E-04*(G.node[idwall]['length']/2*G.node[idwall]['height'])/G.node[idwall]['thickness']
                    matrix[idwall][idwall] -= K #Doussan matrix
                    rhs_p[idwall][1] = -K     #Right-hand side vector
                    #rhs_p[idwall] -= Psi_xyl[0][count]*K     #Right-hand side vector
                    listxy.append(idwall)
                    K_xyl.append(K)
        
        #Adding matrix components at junction-xylem connections
        for cle, valeur in dico.items(): #Loop on junctions between walls
            for cle3, valeur3 in dico3.items():
                if cle3 == cle: #Both cle and cle3 contain the X Y coordinates of wall junctions
                    azer=0
                    length=0.0
                    to = valeur3[0] #junction ID
                    for x in valeur: #Loop on walls around the junction. valeur is the list of cell walls ID meeting at the junction cle
                        if int(x) in listxy: #if the wall around the junction is at interface with xylem
                            azer +=1
                            length+=lengths[int(x)]/4.0 #If it's a border junction wall, 25% of each neighbouring border cell wall is attributed to the length of this jucntion wall
                    if azer>=2: #if the junction connects two walls at interface with xylem
                        K=kw*1.0E-04*length*G.node[to]['height']/G.node[to]['thickness']
                        matrix[to][to] -= K #Doussan matrix
                        rhs_p[to][1] = -K #Right-hand side vector
                        #rhs_p[to] -= Psi_xyl[0][count]*K #Right-hand side vector
                        listxy.append(to)
                        K_xyl.append(K)
    				
        #Adding matrix components at soil-junction connections
        for cle, valeur in dico.items(): #Loop on junctions between walls
            for cle3, valeur3 in dico3.items():
                if cle3 == cle: #Both cle and cle3 contain the X Y coordinates of wall junctions
                    azer=0
                    length=0.0
                    to = valeur3[0] #junction ID
                    for x in valeur: #Loop on walls around the junction. valeur is the list of cell walls ID meeting at the junction cle
                        if int(x) in listesol: #if the wall around the junction is at interface with soil
                            azer +=1
                            length+=lengths[int(x)]/4.0 #If it's a border junction wall, 25% of each neighbouring border cell wall is attributed to the length of this jucntion wall
                    if azer>=2: #if the junction connects two walls at interface with soil
                        K=kw*1.0E-04*length*G.node[to]['height']/G.node[to]['thickness']
                        matrix[to][to] -= K #Doussan matrix
                        rhs_p[to][0] = -K #Right-hand side vector
                        #rhs_p[to] -= Psi_soil[0][count]*K#Right-hand side vector
                        listesol.append(to)
                        K_soil.append(K)
        
        rhs = np.dot(rhs_p,array([[Psi_soil[0][count]],[Psi_xyl[0][count]]])) + rhs_o #multiplication of rhs_p delayed till this point so that rhs_p can be re-used
        #rhs = rhs_p + rhs_o
        ##################################################
        ##Solve Doussan equation, results in soln matrix##
        ##################################################
    
        soln = np.linalg.solve(matrix,rhs) #Solving the equation to get potentials inside the network
    
        #Verification that computation was correct
        verif1=np.allclose(np.dot(matrix,soln),rhs)
    
        #print("Correct computation on PSI ?", verif1)
    
        #Flow rates at interfaces
        Q_soil=[]
        i=0
        for ind in listesol:
            Q_soil.append(K_soil[i]*(Psi_soil[0][count]-soln[ind])) #(cm^3/d) Positive for water flowing into the root
            i+=1
        Q_xyl=[]
        i=0
        for ind in listxy:
            Q_xyl.append(K_xyl[i]*(Psi_xyl[0][count]-soln[ind])) #(cm^3/d) Negative for water flowing into xylem tubes
            cellid=int(Wall2Cell1[ind])
            rank=int(Cell_rank[cellid])
            row=int(rank2row[rank])
            Q_xyl_layer[row][iBarriers][count] += K_xyl[i]*(Psi_xyl[0][count]-soln[ind])
            i+=1
        Q_tot[iBarriers][0]=sum(Q_soil) #(cm^3/d)
        kr_tot[iBarriers][0]=Q_tot[iBarriers][0]/(Psi_soil[0][0]-Psi_xyl[0][0])/perimeter/height/1.0E-04
        print("S0, flow rate per unit root length: soil ",(sum(Q_soil)/height/1.0E-04),"cm^2/d, xylem ",(sum(Q_xyl)/height/1.0E-04),"cm^2/d")
        print("Mass balance error:",(sum(Q_soil)+sum(Q_xyl))/height/1.0E-04,"cm^2/d")
        print("Radial conductivity:",kr_tot[iBarriers][0],"cm/hPa/d, barrier:",Maturity)
        
        #Calculation of standard transmembrane fractions
        jmb=0 #Index for membrane conductance vector
        for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
            i = indice[node] #Node ID number
            psi = soln[i]    #Node water potential
            psi_o_cell = inf #Opposite cell water potential
            count_endo=0 #total number of endodermis cells around the wall
            count_peri=0 #total number of pericycle cells around the wall
            count_exo=0 #total number of exodermis cells around the wall
            count_epi=0 #total number of epidermis cells around the wall
            #count_stele=0 #total number of epidermis cells around the wall
            count_cortex=0 #total number of epidermis cells around the wall
            count_passage=0 #total number of passage cells around the wall
            if i<coro: #wall ID
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == 'membrane': #Wall connection
                        if any(passage_cell_ID==array((indice[neighboor])-zoro)):
                            count_passage+=1
                        if G.node[neighboor]['cgroup']==3:#Endodermis
                            count_endo+=1
                        elif G.node[neighboor]['cgroup']==16 or G.node[neighboor]['cgroup']==5:#Pericycle or stele
                            count_peri+=1
                        elif G.node[neighboor]['cgroup']==1:#Exodermis
                            count_exo+=1
                        elif G.node[neighboor]['cgroup']==2:#Epidermis
                            count_epi+=1
                        elif G.node[neighboor]['cgroup']==4:#Cortex
                            count_cortex+=1
                    # if G.node[neighboor]['cgroup']==5:#Stele
                    #     count_stele+=1
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    j = indice[neighboor] #Neighbouring node ID number
                    psin = soln[j] #Neighbouring node water potential
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if path == "membrane": #Membrane connection
                        K=Kmb[jmb]
                        jmb+=1
                        #Flow densities calculation
                        #Macroscopic distributed parameter for transmembrane flow
                        #Discretization based on cell layers and apoplasmic barriers
                        rank = int(Cell_rank[j-zoro])
                        row = int(rank2row[rank])
                        if rank == 1 and count_epi > 0: #Outer exodermis
                            row += 1
                        if rank == 3 and count_cortex > 0: #Outer endodermis
                            if any(passage_cell_ID==array(j-zoro)) and Maturity==2:
                                row += 2
                            else:
                                row += 3
                        elif rank == 3 and count_peri > 0: #Inner endodermis
                            if any(passage_cell_ID==array(j-zoro)) and Maturity==2:
                                row += 1
                        Flow = K * (psi - psin)
                        if j-zoro not in InterCid: #Not part of STF if crosses an intercellular space "membrane" (that is no membrane though still labelled like one)
                            if Flow > 0 :
                                UptakeLayer_plus[row][iBarriers][count] += Flow #grouping membrane flow rates in cell layers
                            else:
                                UptakeLayer_minus[row][iBarriers][count] += Flow
                            if Flow/Q_tot[iBarriers][0] > 0 :
                                STFlayer_plus[row][iBarriers] += Flow/Q_tot[iBarriers][0] #Cell standard transmembrane fraction (positive)
                                STFcell_plus[j-zoro][iBarriers] += Flow/Q_tot[iBarriers][0] #Cell standard transmembrane fraction (positive)
                            else:
                                STFlayer_minus[row][iBarriers] += Flow/Q_tot[iBarriers][0] #Cell standard transmembrane fraction (negative)
                                STFcell_minus[j-zoro][iBarriers] += Flow/Q_tot[iBarriers][0] #Cell standard transmembrane fraction (negative)
        
        for count in range(1,Nscenarios):
            
            #Soil and xylem water potentials
            Psi_soil[0][count]=float(Psi_soil_range[count].get("pressure")) #Soil water potential (hPa)
            Psi_xyl[0][count]=float(Psi_xyl_range[count].get("pressure")) #Soil water potential (hPa)
            
            ##Reflection coefficients of membranes (undimensional)
            s_factor[0][count]=float(Psi_cell_range[count].get("s_factor")) #(undimensional [0 -> 1]) multiplies all sigma values
            s_membrane=s_factor[0][count]*1.0
            
            #Osmotic potentials (hPa)
            Os_cortex[0][count]=float(Psi_cell_range[count].get("Os_cortex")) # Cortical cell osmotic potential (hPa)
            Os_soil[0][count]=float(Psi_soil_range[count].get("osmotic"))
            Os_xyl[0][count]=float(Psi_xyl_range[count].get("osmotic"))
            #Os_apo=-3000 #-0.3 MPa (Enns et al., 2000) applied stress
            #-0.80 MPa (Enns et al., 2000) concentration of cortical cells, no KNO3
            Os_cell=float(Os_cortex[0][count])
            
            rhs_o = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of osmotic potentials
            jmb=0 #Index for membrane conductance vector
            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                i=indice[node] #Node ID number
                count_endo=0 #total number of endodermis cells around the wall
                count_peri=0 #total number of pericycle cells around the wall
                count_exo=0 #total number of exodermis cells around the wall
                count_epi=0 #total number of epidermis cells around the wall
                count_cortex=0 #total number of cortical cells around the wall
                count_passage=0 #total number of passage cells around the wall
                if i<coro: #wall ID
                    Os_soil_local=float(Os_soil[0][count])
                    for neighboor, eattr in edges.items(): #Loop on connections (edges)
                        if eattr['path'] == 'membrane': #Wall connection
                            if any(passage_cell_ID==array((indice[neighboor])-zoro)):
                                count_passage+=1
                            if G.node[neighboor]['cgroup']==3:#Endodermis
                                count_endo+=1
                            elif G.node[neighboor]['cgroup']==16 or G.node[neighboor]['cgroup']==5:#Pericycle or stele
                                count_peri+=1
                            elif G.node[neighboor]['cgroup']==4:#Cortex
                                count_cortex+=1
                            elif G.node[neighboor]['cgroup']==1:#Exodermis
                                count_exo+=1
                            elif G.node[neighboor]['cgroup']==2:#Epidermis
                                count_epi+=1
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    j = (indice[neighboor]) #neighbouring node number
                    if j > i: #Only treating the information one way to save time
                        path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                        if path == "membrane": #Membrane connection
                            #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                            if G.node[j]['cgroup']==1:#Exodermis
                                Os_wall=Os_soil_local
                            elif G.node[j]['cgroup']==2:#Epidermis
                                Os_wall=Os_soil_local
                            elif G.node[j]['cgroup']==3:#Endodermis
                                if count_peri==0: #wall between endodermis and cortex or between two endodermal cells
                                    Os_wall=Os_soil_local
                                elif count_peri>0 and count_endo>0: #wall between endodermis and pericycle
                                    if Maturity==0: #No apoplastic barrier
                                        Os_wall=Os_soil_local
                                    else:
                                        Os_wall=float(Os_xyl[0][count])
                            elif G.node[j]['cgroup']==4:#Cortex
                                Os_wall=Os_soil_local
                            elif G.node[j]['cgroup']==5:#Stele
                                if Maturity==0: #No apoplastic barrier
                                    Os_wall=Os_soil_local
                                else:
                                    Os_wall=float(Os_xyl[0][count])
                            elif G.node[j]['cgroup']==16:#Pericycle
                                if Maturity==0: #No apoplastic barrier
                                    Os_wall=Os_soil_local
                                else:
                                    Os_wall=float(Os_xyl[0][count])
                            K=Kmb[jmb]
                            jmb+=1
                            rhs_o[i]+= K*s_membrane*(Os_wall - Os_cell) #Wall node
                            rhs_o[j]+= K*s_membrane*(Os_cell - Os_wall) #Cell node
    
            rhs = np.dot(rhs_p,array([[Psi_soil[0][count]],[Psi_xyl[0][count]]])) + rhs_o #multiplication of rhs_p delayed till this point so that rhs_p can be re-used
            #rhs = rhs_p + rhs_o
            ##################################################
            ##Solve Doussan equation, results in soln matrix##
            ##################################################
            
            #soln = slin.solve(matrix,rhs)
            soln = np.linalg.solve(matrix,rhs) #Solving the equation to get potentials inside the network
            
            #Verification that computation was correct
            verif1=np.allclose(np.dot(matrix,soln),rhs)
            
            #print("Correct computation on PSI ?", verif1)
            
            #Flow rates at interfaces
            Q_soil=[]
            i=0
            for ind in listesol:
                Q_soil.append(K_soil[i]*(Psi_soil[0][count]-soln[ind])) #(cm^3/d) Positive for water flowing into the root
                i+=1
            Q_xyl=[]
            i=0
            for ind in listxy:
                Q_xyl.append(K_xyl[i]*(Psi_xyl[0][count]-soln[ind])) #(cm^3/d) Negative for water flowing into xylem tubes
                cellid=int(Wall2Cell1[ind])
                rank=int(Cell_rank[cellid])
                row=int(rank2row[rank])
                Q_xyl_layer[row][iBarriers][count] += K_xyl[i]*(Psi_xyl[0][count]-soln[ind])
                i+=1
            Q_tot[iBarriers][count]=sum(Q_soil)
            for ind in range(zoro,len(G.node)): #zoro is the index of the first cell
                cellnumber1=ind-zoro
                rank = int(Cell_rank[cellnumber1])
                row = int(rank2row[rank])
                if rank == 1: #Exodermis
                    PsiCellLayer[row][iBarriers][count] += soln[ind]*(STFcell_plus[cellnumber1][iBarriers]+abs(STFcell_minus[cellnumber1][iBarriers]))/(STFlayer_plus[row][iBarriers]+abs(STFlayer_minus[row][iBarriers])+STFlayer_plus[row+1][iBarriers]+abs(STFlayer_minus[row+1][iBarriers])) #(hPa)
                    PsiCellLayer[row+1][iBarriers][count] += soln[ind]*(STFcell_plus[cellnumber1][iBarriers]+abs(STFcell_minus[cellnumber1][iBarriers]))/(STFlayer_plus[row][iBarriers]+abs(STFlayer_minus[row][iBarriers])+STFlayer_plus[row+1][iBarriers]+abs(STFlayer_minus[row+1][iBarriers])) #(hPa)
                elif rank == 3: #Endodermis
                    if any(passage_cell_ID==array(cellnumber1)) and Maturity==2: #Passage cell
                        PsiCellLayer[row+1][iBarriers][count] += soln[ind]*(STFcell_plus[cellnumber1][iBarriers]+abs(STFcell_minus[cellnumber1][iBarriers]))/(STFlayer_plus[row+1][iBarriers]+abs(STFlayer_minus[row+1][iBarriers])+STFlayer_plus[row+2][iBarriers]+abs(STFlayer_minus[row+2][iBarriers])) #(hPa)
                        PsiCellLayer[row+2][iBarriers][count] += soln[ind]*(STFcell_plus[cellnumber1][iBarriers]+abs(STFcell_minus[cellnumber1][iBarriers]))/(STFlayer_plus[row+1][iBarriers]+abs(STFlayer_minus[row+1][iBarriers])+STFlayer_plus[row+2][iBarriers]+abs(STFlayer_minus[row+2][iBarriers])) #(hPa)
                    else:
                        PsiCellLayer[row][iBarriers][count] += soln[ind]*(STFcell_plus[cellnumber1][iBarriers]+abs(STFcell_minus[cellnumber1][iBarriers]))/(STFlayer_plus[row][iBarriers]+abs(STFlayer_minus[row][iBarriers])+STFlayer_plus[row+3][iBarriers]+abs(STFlayer_minus[row+3][iBarriers])) #(hPa)
                        PsiCellLayer[row+3][iBarriers][count] += soln[ind]*(STFcell_plus[cellnumber1][iBarriers]+abs(STFcell_minus[cellnumber1][iBarriers]))/(STFlayer_plus[row][iBarriers]+abs(STFlayer_minus[row][iBarriers])+STFlayer_plus[row+3][iBarriers]+abs(STFlayer_minus[row+3][iBarriers])) #(hPa)
                        if not Maturity==2:
                            PsiCellLayer[row+1][iBarriers][count] = nan
                            PsiCellLayer[row+2][iBarriers][count] = nan
                else:
                    PsiCellLayer[row][iBarriers][count] += soln[ind]*(STFcell_plus[cellnumber1][iBarriers]+abs(STFcell_minus[cellnumber1][iBarriers]))/(STFlayer_plus[row][iBarriers]+abs(STFlayer_minus[row][iBarriers])) #(hPa)
            i=0
            for ind in rank2row:
                if not isnan(ind):
                    ind=int(ind)
                    if i==2: #Epidermis
                        OsCellLayer[ind][count]=Os_cell
                        OsWallLayer[ind][count]=Os_soil[0][count]  #Not valid any more
                    elif i==1: #Exodermis
                        OsCellLayer[ind][count]=Os_cell
                        OsCellLayer[ind+1][count]=Os_cell
                        OsWallLayer[ind][count]=Os_soil[0][count]
                        OsWallLayer[ind+1][count]=Os_soil[0][count]
                    elif i==3: #Endodermis
                        OsCellLayer[ind][count]=Os_cell
                        OsCellLayer[ind+3][count]=Os_cell
                        if Maturity==2:
                            OsCellLayer[ind+1][count]=Os_cell
                            OsCellLayer[ind+2][count]=Os_cell
                        else:
                            OsCellLayer[ind+1][count]=nan
                            OsCellLayer[ind+2][count]=nan
                        OsWallLayer[ind][count]=Os_xyl[0][count]
                        OsWallLayer[ind+3][count]=Os_soil[0][count]
                        if Maturity==2:
                            OsWallLayer[ind+1][count]=Os_xyl[0][count]
                            OsWallLayer[ind+2][count]=Os_soil[0][count]
                        else:
                            OsWallLayer[ind+1][count]=nan
                            OsWallLayer[ind+2][count]=nan
                    elif i==16: #Pericycle
                        OsCellLayer[ind][count]=Os_cell
                        OsWallLayer[ind][count]=Os_xyl[0][count]
                    elif i>=40 and i<50: #Cortex
                        OsCellLayer[ind][count]=Os_cell
                        OsWallLayer[ind][count]=Os_soil[0][count]
                    elif i>=50: #Stele
                        OsCellLayer[ind][count]=Os_cell
                        OsWallLayer[ind][count]=Os_xyl[0][count]
                i+=1
            print("Flow rate per unit root length: soil ",(sum(Q_soil)/height/1.0E-04),"cm^2/d, xylem ",(sum(Q_xyl)/height/1.0E-04),"cm^2/d")
            print("Mass balance error:",(sum(Q_soil)+sum(Q_xyl))/height/1.0E-04,"cm^2/d")
            
            #################################################################
            ##Calul of Fluxes between nodes and Creating the edge_flux_list##
            #################################################################
            
            #Creating a list for the fluxes
            #edge_flux_list=[]
            
            #Filling the fluxes list
            MembraneFlowDensity=[]
            WallFlowDensity=[]
            PlasmodesmFlowDensity=[]
            jmb=0 #Index for membrane conductance vector
            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                i = indice[node] #Node ID number
                psi = soln[i] #Node water potential
                psi_o_cell = inf #Opposite cell water potential
                count_endo=0 #total number of endodermis cells around the wall
                count_peri=0 #total number of pericycle cells around the wall
                count_exo=0 #total number of exodermis cells around the wall
                count_epi=0 #total number of epidermis cells around the wall
                count_stele=0 #total number of epidermis cells around the wall
                count_cortex=0 #total number of epidermis cells around the wall
                count_passage=0 #total number of passage cells around the wall
                noPD=False #Initializes the flag for wall connected to an intercellular space -> does not have plasmodesmata
                if i<coro: #wall ID
                    Os_soil_local=float(Os_soil[0][count])
                    for neighboor, eattr in edges.items(): #Loop on connections (edges)
                        if eattr['path'] == 'membrane': #Wall connection
                            if any(passage_cell_ID==array((indice[neighboor])-zoro)):
                                count_passage+=1
                            if G.node[neighboor]['cgroup']==3:#Endodermis
                                count_endo+=1
                            elif G.node[neighboor]['cgroup']==16 or G.node[neighboor]['cgroup']==5:#Pericycle or stele
                                count_peri+=1
                            elif G.node[neighboor]['cgroup']==1:#Exodermis
                                count_exo+=1
                            elif G.node[neighboor]['cgroup']==2:#Epidermis
                                count_epi+=1
                            elif G.node[neighboor]['cgroup']==4:#Cortex
                                count_cortex+=1
                            if G.node[neighboor]['cgroup']==5:#Stele
                                count_stele+=1
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    j = indice[neighboor] #Neighbouring node ID number
                    #if j > i: #Only treating the information one way to save time
                    psin = soln[j] #Neighbouring node water potential
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if i<coro:
                        if Paraview==1:
                            if path == "wall":
                                #K = eattr['kw']*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                if count_endo>=2: #wall between two endodermis cells
                                    K = kw_endo_endo*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                elif count_peri>0 and count_endo>0: #wall between endodermis and pericycle
                                    if count_passage>0:
                                        K = kw_passage*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length']
                                    else:
                                        K = kw_endo_peri*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                elif count_peri==0 and count_endo==1: #wall between endodermis and pericycle
                                    if count_passage>0:
                                        K = kw_passage*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length']
                                    else:
                                        K = kw_endo_cortex*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                elif count_exo>=2: #wall between two exodermis cells
                                    K = kw_exo_exo*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                else: #other walls
                                    K = kw*1.0E-04*((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                #The ordering in WallFlowDensity will correspond to the one of ThickWallsX, saved for display only
                                WallFlowDensity.append((i,j,K / (((eattr['lat_dist']+eattr['height'])*eattr['thickness']-square(eattr['thickness']))*1.0E-08) * (psin - psi) * sign(j-i))) # (cm/d) Positive towards lower node ID 
                                #WallFlowDensity.append((i,j,eattr['kw']/(eattr['length']*1.0E-04) * (psin - psi) * sign(j-i))) #Positive towards lower node ID 
                            elif path == "membrane": #Membrane connection
                                #K = (eattr['kmb']+eattr['kaqp'])*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                                if j-zoro in InterCid: #the neighbour is an intercellular space "cell"
                                    kaqp=kInterC
                                    noPD=True
                                else:
                                    kaqp=kaqp_cortex
                                if count_endo>=2: #wall between two endodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if kw_endo_endo==0.00:
                                        K=0.00
                                    else:
                                        K = 1/(1/(kw_endo_endo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length'] #(cm^3/hPa/d)
                                elif count_exo>=2: #wall between two exodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if kw_exo_exo==0.00:
                                        K=0.00
                                    else:
                                        K = 1/(1/(kw_exo_exo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length'] #(cm^3/hPa/d)
                                elif count_peri>0 and count_endo>0: #wall between endodermis and pericycle, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if count_passage>0:
                                        K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                                    else:
                                        if kw_endo_peri==0.00:
                                            K=0.00
                                        else:
                                            K = 1/(1/(kw_endo_peri/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                                elif count_peri==0 and count_endo==1: #wall between endodermis and pericycle, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if kaqp==0.0:
                                        K=1.00E-14
                                    else:
                                        if count_passage>0:
                                            K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                                        else:
                                            if kw_endo_cortex==0.00:
                                                K=0.00
                                            else:
                                                K = 1/(1/(kw_endo_cortex/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']
                                else:
                                    if kaqp==0.0:
                                        K=1.00E-14
                                    else:
                                        K = 1/(1/(kw/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(eattr['height']+eattr['dist'])*eattr['length'] #(cm^3/hPa/d)
                                #Flow densities calculation
                                #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                                if G.node[j]['cgroup']==1:#Exodermis
                                    Os_wall=Os_soil_local
                                elif G.node[j]['cgroup']==2:#Epidermis
                                    Os_wall=Os_soil_local
                                elif G.node[j]['cgroup']==3:#Endodermis
                                    if count_peri==0: #wall between endodermis and cortex or between two endodermal cells
                                        Os_wall=Os_soil_local
                                    elif count_peri>0: #wall between endodermis and pericycle
                                        if Maturity==0: #No apoplastic barrier
                                            Os_wall=Os_soil_local
                                        else:
                                            Os_wall=float(Os_xyl[0][count])
                                elif G.node[j]['cgroup']==4:#Cortex
                                    Os_wall=Os_soil_local
                                elif G.node[j]['cgroup']==5:#Stele
                                    if Maturity==0: #No apoplastic barrier
                                        Os_wall=Os_soil_local
                                    else:
                                        Os_wall=float(Os_xyl[0][count])
                                elif G.node[j]['cgroup']==16:#Pericycle
                                    if Maturity==0: #No apoplastic barrier
                                        Os_wall=Os_soil_local
                                    else:
                                        Os_wall=float(Os_xyl[0][count])
                                #The ordering in MembraneFlowDensity will correspond to the one of ThickWalls, saved for display only
                                MembraneFlowDensity.append(K / (1.0E-08*(eattr['height']+eattr['dist'])*eattr['length']) * (psi - psin + s_membrane*(Os_wall - Os_cell)))
                                #Macroscopic distributed parameter for transmembrane flow
                                #Discretization based on cell layers and apoplasmic barriers
                                rank = int(Cell_rank[j-zoro])
                                row = int(rank2row[rank])
                                if rank == 1 and count_epi > 0: #Outer exodermis
                                    row += 1
                                if rank == 3 and count_cortex > 0: #Outer endodermis
                                    if any(passage_cell_ID==array(j-zoro)) and Maturity==2:
                                        row += 2
                                    else:
                                        row += 3
                                elif rank == 3 and count_peri > 0: #Inner endodermis
                                    if any(passage_cell_ID==array(j-zoro)) and Maturity==2:
                                        row += 1
                                        #print('PsiWallPassage:',psi)
                                Flow = K * (psi - psin + s_membrane*(Os_wall - Os_cell))
                                if j-zoro not in InterCid:
                                    if Flow > 0 :
                                        UptakeLayer_plus[row][iBarriers][count] += Flow #grouping membrane flow rates in cell layers
                                    else:
                                        UptakeLayer_minus[row][iBarriers][count] += Flow
                                
                                if K>1.0e-12: #Not an impermeable wall
                                    PsiWallLayer[row][iBarriers][count] += psi
                                    NWallLayer[row][iBarriers][count] += 1
                                
                                if psi_o_cell == inf:
                                    psi_o_cell=psin
                                else:
                                    if noPD: #No plasmodesmata because the wall i is connected to an intercellular space
                                        PlasmodesmFlowDensity.append(0.0)
                                        PlasmodesmFlowDensity.append(0.0) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif count_epi==1 and count_exo==1: #wall between epidermis and exodermis
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_epi_exo * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_epi_exo * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif (count_exo==1 or count_epi==1) and count_cortex==1: #wall between exodermis and cortex
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_outer_cortex * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_outer_cortex * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif count_cortex==2: #wall between cortical cells
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_cortex_cortex * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_cortex_cortex * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif count_cortex==1 and count_endo==1: #wall between cortex and endodermis
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_peri_stele * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_peri_stele * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif count_endo==2: #wall between endodermal cells
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_endo_endo * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_endo_endo * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif count_peri>0 and count_endo>0: #wall between endodermis and pericycle
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_endo_peri * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_endo_peri * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif count_stele==2: #wall between stele cells
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_stele_stele * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_stele_stele * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    elif count_peri>0 and count_stele==1: #wall between stele and pericycle
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_peri_stele * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight_peri_stele * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                    else: #Default plasmodesmatal frequency
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight * (psin - psi_o_cell) / (1.0E-04*eattr['height']))
                                        PlasmodesmFlowDensity.append(Kpl*Fplxheight * (psi_o_cell - psin) / (1.0E-04*eattr['height'])) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only 
                        elif Paraview==0:
                            if path == "membrane": #Membrane connection
                                K=Kmb[jmb]
                                jmb+=1
                                #Flow densities calculation
                                #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                                if G.node[j]['cgroup']==1:#Exodermis
                                    Os_wall=Os_soil_local
                                elif G.node[j]['cgroup']==2:#Epidermis
                                    Os_wall=Os_soil_local
                                elif G.node[j]['cgroup']==3:#Endodermis
                                    if count_peri==0: #wall between endodermis and cortex or between two endodermal cells
                                        Os_wall=Os_soil_local
                                    elif count_peri>0: #wall between endodermis and pericycle
                                        if Maturity==0: #No apoplastic barrier
                                            Os_wall=Os_soil_local
                                        else:
                                            Os_wall=float(Os_xyl[0][count])
                                elif G.node[j]['cgroup']==4:#Cortex
                                    Os_wall=Os_soil_local
                                elif G.node[j]['cgroup']==5:#Stele
                                    if Maturity==0: #No apoplastic barrier
                                        Os_wall=Os_soil_local
                                    else:
                                        Os_wall=float(Os_xyl[0][count])
                                elif G.node[j]['cgroup']==16:#Pericycle
                                    if Maturity==0: #No apoplastic barrier
                                        Os_wall=Os_soil_local
                                    else:
                                        Os_wall=float(Os_xyl[0][count])
                                #Macroscopic distributed parameter for transmembrane flow
                                #Discretization based on cell layers and apoplasmic barriers
                                rank = int(Cell_rank[j-zoro])
                                row = int(rank2row[rank])
                                if rank == 1 and count_epi > 0: #Outer exodermis
                                    row += 1
                                if rank == 3 and count_cortex > 0: #Outer endodermis
                                    if any(passage_cell_ID==array(j-zoro)) and Maturity==2:
                                        row += 2
                                    else:
                                        row += 3
                                elif rank == 3 and count_peri > 0: #Inner endodermis
                                    if any(passage_cell_ID==array(j-zoro)) and Maturity==2:
                                        row += 1
                                Flow = K * (psi - psin + s_membrane*(Os_wall - Os_cell))
                                if j-zoro not in InterCid:
                                    if Flow > 0 :
                                        UptakeLayer_plus[row][iBarriers][count] += Flow #grouping membrane flow rates in cell layers
                                    else:
                                        UptakeLayer_minus[row][iBarriers][count] += Flow
                                
                                if K>1.0e-12: #Not an impermeable wall
                                    PsiWallLayer[row][iBarriers][count] += psi
                                    NWallLayer[row][iBarriers][count] += 1
            
            ####################################
            ## Creates .vtk file for Paraview ##
            ####################################
            
            if Paraview==1:
                if ParaviewWP==1: #2D visualization of walls pressure potentials
                    text_file = open(newpath+"Walls2Db"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    sath0=max(soln[0:zoro-1])
                    satl0=min(soln[0:zoro-1])
                    with open(newpath+"Walls2Db"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Wall geometry 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(coro*2) + " " + str(coro*6) + " \n") #len(G.node)
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                j=indice[neighboor]
                                if j>i and eattr['path']=='wall':
                                    #print(nx.get_node_attributes(edges,'path'))
                                    myfile.write(str(2) + " " + str(i) + " " + str(j) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(coro*2) + " \n") #The number of nodes corresponds to the number of wall to wall connections.... to be checked, might not be generality
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                j=indice[neighboor]
                                if j>i and eattr['path']=='wall':
                                    #print(nx.get_node_attributes(edges,'path'))
                                    myfile.write(str(3) + " \n") #Line cell type
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                        myfile.write("SCALARS Wall_pressure float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for node in G:
                            myfile.write(str(float(min(sath0,max(satl0,soln[node])))) + " \n") #Line cell type
                    myfile.close()
                    text_file.close()
                
                if ParaviewWP==1: #2D visualization of walls osmotic potentials
                    text_file = open(newpath+"WallsOs2Db"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"WallsOs2Db"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Wall geometry 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(coro*2) + " " + str(coro*6) + " \n") #len(G.node)
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                j=indice[neighboor]
                                if j>i and eattr['path']=='wall':
                                    #print(nx.get_node_attributes(edges,'path'))
                                    myfile.write(str(2) + " " + str(i) + " " + str(j) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(coro*2) + " \n") #The number of nodes corresponds to the number of wall to wall connections.... to be checked, might not be generality
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                j=indice[neighboor]
                                if j>i and eattr['path']=='wall':
                                    #print(nx.get_node_attributes(edges,'path'))
                                    myfile.write(str(3) + " \n") #Line cell type
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                        myfile.write("SCALARS Wall_osmotic_pot float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for node, edges in G.adjacency_iter():
                            i=indice[node] #Node ID number
                            count_endo=0 #total number of endodermis cells around the wall
                            count_peri=0 #total number of pericycle cells around the wall
                            count_exo=0 #total number of exodermis cells around the wall
                            count_epi=0 #total number of epidermis cells around the wall
                            count_cortex=0 #total number of cortical cells around the wall
                            count_passage=0 #total number of passage cells around the wall
                            if i<coro: #wall ID
                                Os_soil_local=float(Os_soil[0][count])
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    if eattr['path'] == 'membrane': #Wall connection
                                        if G.node[neighboor]['cgroup']==3:#Endodermis
                                            count_endo+=1
                                        elif G.node[neighboor]['cgroup']==16 or G.node[neighboor]['cgroup']==5:#Pericycle or stele
                                            count_peri+=1
                                        elif G.node[neighboor]['cgroup']==4:#Cortex
                                            count_cortex+=1
                                        elif G.node[neighboor]['cgroup']==1:#Exodermis
                                            count_exo+=1
                                        elif G.node[neighboor]['cgroup']==2:#Epidermis
                                            count_epi+=1
                            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                j = (indice[neighboor]) #neighbouring node number
                                if j > i: #Only treating the information one way to save time
                                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                                    if path == "membrane": #Membrane connection
                                        #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                                        if G.node[j]['cgroup']==1:#Exodermis
                                            Os_wall=Os_soil_local
                                        elif G.node[j]['cgroup']==2:#Epidermis
                                            Os_wall=Os_soil_local
                                        elif G.node[j]['cgroup']==3:#Endodermis
                                            if count_peri==0: #wall between endodermis and cortex or between two endodermal cells
                                                Os_wall=Os_soil_local
                                            elif count_peri>0: #wall between endodermis and pericycle
                                                if Maturity==0: #No apoplastic barrier
                                                    Os_wall=Os_soil_local
                                                else:
                                                    Os_wall=float(Os_xyl[0][count])
                                        elif G.node[j]['cgroup']==4:#Cortex
                                            Os_wall=Os_soil_local
                                        elif G.node[j]['cgroup']==5:#Stele
                                            if Maturity==0: #No apoplastic barrier
                                                Os_wall=Os_soil_local
                                            else:
                                                Os_wall=float(Os_xyl[0][count])
                                        elif G.node[j]['cgroup']==16:#Pericycle
                                            if Maturity==0: #No apoplastic barrier
                                                Os_wall=Os_soil_local
                                            else:
                                                Os_wall=float(Os_xyl[0][count])
                            myfile.write(str(float(Os_wall)) + " \n") #Line cell type
                    myfile.close()
                    text_file.close()
            
                
                if ParaviewWP==1 and ParaviewCP==1: #2D visualization of walls & cells water potentials
                    text_file = open(newpath+"WallsAndCells2Db"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"WallsAndCells2Db"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Water potential distribution in cells and walls 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(coro*2+Ncells) + " " + str(coro*6+Ncells*2) + " \n") #
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                j=indice[neighboor]
                                if j>i and eattr['path']=='wall':
                                    #print(nx.get_node_attributes(edges,'path'))
                                    myfile.write(str(2) + " " + str(i) + " " + str(j) + " \n")
                            if i>=zoro: #Cell node
                                myfile.write("1 " + str(i) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(coro*2+Ncells) + " \n") #
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                j=indice[neighboor]
                                if j>i and eattr['path']=='wall':
                                    #print(nx.get_node_attributes(edges,'path'))
                                    myfile.write(str(3) + " \n") #Line cell type
                            if i>=zoro: #Cell node
                                myfile.write("1 \n")
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                        myfile.write("SCALARS pressure float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for node in G:
                            myfile.write(str(float(soln[node])) + " \n") #Line cell type
                    myfile.close()
                    text_file.close()
                
                if ParaviewCP==1: #2D visualization of cells water potentials
                    text_file = open(newpath+"Cells2Db"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    sath01=max(soln[zoro:zoro+Ncells-1])
                    satl01=min(soln[zoro:zoro+Ncells-1])
                    with open(newpath+"Cells2Db"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Pressure potential distribution in cells 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(Ncells) + " " + str(Ncells*2) + " \n") #
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            if i>=zoro: #Cell node
                                myfile.write("1 " + str(i) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(Ncells) + " \n") #
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            if i>=zoro: #Cell node
                                myfile.write("1 \n")
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                        myfile.write("SCALARS Cell_pressure float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for node in G:
                            myfile.write(str(float(min(sath01,max(satl01,soln[node])))) + " \n") #Line cell type
                    myfile.close()
                    text_file.close()
                
                if ParaviewMF==1: #3D visualization of membrane fluxes
                    text_file = open(newpath+"Membranes3Db"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    sath1=max(MembraneFlowDensity)*color_threshold
                    satl1=min(MembraneFlowDensity)*color_threshold
                    if satl1<-sath1: #min(MembraneFlowDensity)<0:
                        sath1=-satl1
                    else:
                        satl1=-sath1
                    with open(newpath+"Membranes3Db"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Membranes geometry 3D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWalls)*2)+" float \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " 0.0 \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(len(ThickWalls)) + " " + str(len(ThickWalls)*5) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[1]>=coro: #wall that is a junction
                                myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[5])) + " " + str(int(ThickWallNode[5])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n") #All points were repeated twice (once at z=0 and once at z=height), so adding len(ThickWalls) is the same point at z=height
                                myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[6])) + " " + str(int(ThickWallNode[6])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(len(ThickWalls)) + " \n")
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[1]>=coro: #wall that is a junction
                                myfile.write("9 \n") #Quad cell type
                                myfile.write("9 \n") #Quad cell type
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(ThickWalls)*2) + " \n")
                        myfile.write("SCALARS TM_flux_(m/s) float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[0]<len(MembraneFlowDensity):
                                myfile.write(str(float(min(sath1,max(satl1,MembraneFlowDensity[int(ThickWallNode[0])])))/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                            else:
                                myfile.write(str(float(min(sath1,max(satl1,(MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)))/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[0]<len(MembraneFlowDensity):
                                myfile.write(str(float(min(sath1,max(satl1,MembraneFlowDensity[int(ThickWallNode[0])])))/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                            else:
                                myfile.write(str(float(min(sath1,max(satl1,(MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)))/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                    myfile.close()
                    text_file.close()
                
                if ParaviewWF==1: #Wall flow density data
                    maxWallFlowDensity=0.0
                    for ir in range(len(WallFlowDensity)):
                        maxWallFlowDensity=max(maxWallFlowDensity,abs(WallFlowDensity[ir][2]))
                    sath2=maxWallFlowDensity*color_threshold #(1-(1-color_threshold)/2)
                    #satl2=0.0
                    text_file = open(newpath+"WallsThick3D_bottomb"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"WallsThick3D_bottomb"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Wall geometry 3D including thickness bottom \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                        for ThickWallNodeX in ThickWallsX:
                            myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " 0.0 \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(int(cr+coro)) + " " + str(int(2*coro*5+sum(nWall2NewWallX[coro:])+cr-coro+2*len(Wall2NewWallX[coro:]))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                        for PolygonX in ThickWallPolygonX:
                            myfile.write("4 " + str(int(PolygonX[0])) + " " + str(int(PolygonX[1])) + " " + str(int(PolygonX[2])) + " " + str(int(PolygonX[3])) + " \n")
                        i=coro
                        for PolygonX in Wall2NewWallX[coro:]: #"junction" polygons
                            #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                            string=str(int(nWall2NewWallX[i]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                            for id1 in range(nWall2NewWallX[i]):
                                string=string+" "+str(int(PolygonX[id1]))
                            string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                            myfile.write(string + " \n")
                            i+=1
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(cr+coro) + " \n")
                        for PolygonX in ThickWallPolygonX:
                            myfile.write("7 \n") #Polygon cell type
                        for PolygonX in Wall2NewWallX[coro:]:
                            myfile.write("6 \n") #Triangle-strip cell type
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(ThickWallsX)) + " \n")
                        myfile.write("SCALARS Apo_flux_(m/s) float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        NewWallFlowDensity=zeros((len(ThickWallsX),2))
                        i=0
                        for PolygonX in ThickWallPolygonX:
                            for id1 in range(4):
                                if abs(int(WallFlowDensity[i][2]))>min(NewWallFlowDensity[int(PolygonX[id1])]):
                                    NewWallFlowDensity[int(PolygonX[id1])][0]=max(NewWallFlowDensity[int(PolygonX[id1])])
                                    NewWallFlowDensity[int(PolygonX[id1])][1]=abs(float(WallFlowDensity[i][2]))
                            i+=1
                        for i in range(len(ThickWallsX)):
                            #myfile.write(str(min(sath2,max(satl2,log10(mean(NewWallFlowDensity[i]))))) + " \n")
                            myfile.write(str(float(min(sath2,mean(NewWallFlowDensity[i])))/sperd/cmperm) + " \n")
                    myfile.close()
                    text_file.close()
                
                    text_file = open(newpath+"InterC3D_bottomb"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"InterC3D_bottomb"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Intercellular space geometry 3D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height/200) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(len(InterCid)) + " " + str(int(len(InterCid)+sum(nCell2ThickWalls[InterCid]))) + " \n") #The number of cells corresponds to the number of intercellular spaces
                        InterCFlowDensity=zeros((Ncells,1))
                        for cid in InterCid:
                            n=int(nCell2ThickWalls[cid]) #Total number of thick wall nodes around the protoplast
                            Polygon=Cell2ThickWalls[cid][:n]
                            ranking=list()
                            ranking.append(int(Polygon[0]))
                            ranking.append(ThickWalls[int(ranking[0])][5])
                            ranking.append(ThickWalls[int(ranking[0])][6])
                            for id1 in range(1,n):
                                wid1=ThickWalls[int(ranking[id1])][5]
                                wid2=ThickWalls[int(ranking[id1])][6]
                                if wid1 not in ranking:
                                    ranking.append(wid1)
                                if wid2 not in ranking:
                                    ranking.append(wid2)
                            string=str(n)
                            for id1 in ranking:
                                string=string+" "+str(int(id1))
                            myfile.write(string + " \n")
                            for twpid in Polygon[:n/2]: #The first half of nodes are wall nodes actually connected to cells
                                InterCFlowDensity[cellnumber1]+=abs(MembraneFlowDensity[int(twpid)])/n #Mean absolute flow density calculation
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(len(InterCid)) + " \n")
                        for i in range(len(InterCid)):
                            myfile.write("6 \n") #Triangle-strip cell type
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(ThickWalls)) + " \n")
                        myfile.write("SCALARS Apo_flux_(m/s) float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for ThickWallNode in ThickWalls:
                            cellnumber1=ThickWallNode[2]-zoro
                            myfile.write(str(float(min(sath1,max(satl1,InterCFlowDensity[cellnumber1])))/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                    myfile.close()
                    text_file.close()
                
                
                
                if ParaviewPF==1: #Plasmodesmata flow density data disks
                    text_file = open(newpath+"Plasmodesm3Db"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    sath3=max(PlasmodesmFlowDensity)*color_threshold
                    satl3=min(PlasmodesmFlowDensity)*color_threshold
                    if satl3<-sath3: #min(PlasmodesmFlowDensity)<0:
                        sath3=-satl3
                    else:
                        satl3=-sath3
                    with open(newpath+"Plasmodesm3Db"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("PD flux disks 3D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(PlasmodesmFlowDensity)*12)+" float \n")
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[1]<coro: #selection of new walls (not new junctions)
                                if ThickWallNode[7]==0: #new walls that are not at the interface with soil or xylem, where there is no plasmodesmata   #if G.node[int(ThickWallNode[1])]['borderlink']==0
                                    #calculate the XY slope between the two neighbouring new junctions
                                    twpid1=int(ThickWallNode[5])
                                    twpid2=int(ThickWallNode[6])
                                    if not ThickWalls[twpid1][3]==ThickWalls[twpid2][3]: #Otherwise we'll get a division by 0 error
                                        slopeNJ=(ThickWalls[twpid1][4]-ThickWalls[twpid2][4])/(ThickWalls[twpid1][3]-ThickWalls[twpid2][3]) #slope of the line connecting the new junction nodes neighbouring the new wall
                                    else:
                                        slopeNJ=inf
                                    x0=ThickWallNode[3]
                                    y0=ThickWallNode[4]
                                    z0=20.0
                                    #Calculate the horizontal distance between XY0 and the cell center, compare it with the distance between the mean position of the new junctions. If the latter is closer to the cell center, it becomes the new XY0 to make sur the disk is visible
                                    xC=position[int(ThickWallNode[2])][0]
                                    yC=position[int(ThickWallNode[2])][1]
                                    xNJ=(ThickWalls[twpid1][3]+ThickWalls[twpid2][3])/2.0
                                    yNJ=(ThickWalls[twpid1][4]+ThickWalls[twpid2][4])/2.0
                                    if sqrt(square(x0-xC)+square(y0-yC)) > sqrt(square(xNJ-xC)+square(yNJ-yC)):
                                        x0=xNJ
                                        y0=yNJ
                                    for i in range(12):
                                        x=x0+cos(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                        y=y0+sin(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                        z=z0+radiusPlasmodesm_disp*sin(int(i)*pi/6.0)
                                        myfile.write(str(x) + " " + str(y) + " " + str(z) + " \n")
                            else:
                                break #interrupts the for loop in case we reached the new junction nodes
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(len(PlasmodesmFlowDensity)) + " " + str(len(PlasmodesmFlowDensity)*13) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                        for i in range(len(PlasmodesmFlowDensity)):
                            if PlasmodesmFlowDensity[i]==0:
                                myfile.write("12 " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " \n")
                            else:
                                myfile.write("12 " + str(i*12+0) + " " + str(i*12+1) + " " + str(i*12+2) + " " + str(i*12+3) + " " + str(i*12+4) + " " + str(i*12+5) + " " + str(i*12+6) + " " + str(i*12+7) + " " + str(i*12+8) + " " + str(i*12+9) + " " + str(i*12+10) + " " + str(i*12+11) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(len(PlasmodesmFlowDensity)) + " \n")
                        for i in range(len(PlasmodesmFlowDensity)):
                            myfile.write("7 \n") #Polygon cell type 
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(PlasmodesmFlowDensity)*12) + " \n")
                        myfile.write("SCALARS PD_Flux_(m/s) float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for i in range(len(PlasmodesmFlowDensity)):
                            for j in range(12):
                                myfile.write(str(float(min(sath3,max(satl3,PlasmodesmFlowDensity[i])))/sperd/cmperm) + " \n")
                    myfile.close()
                    text_file.close()
                
                
                if ParaviewMF==1 and ParaviewPF==1: #Membranes and plasmodesms in the same file
                    text_file = open(newpath+"Membranes_n_plasmodesm3Db"+str(Maturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"Membranes_n_plasmodesm3Db"+str(Maturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Membranes geometry and plasmodesm disks 3D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWalls)*2+len(PlasmodesmFlowDensity)*12)+" float \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " 0.0 \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height) + " \n")
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[1]<coro: #selection of new walls (not new junctions)
                                if ThickWallNode[7]==0: #new walls that are not at the interface with soil or xylem, where there is no plasmodesmata   #if G.node[int(ThickWallNode[1])]['borderlink']==0
                                    #calculate the XY slope between the two neighbouring new junctions
                                    twpid1=int(ThickWallNode[5])
                                    twpid2=int(ThickWallNode[6])
                                    if not ThickWalls[twpid1][3]==ThickWalls[twpid2][3]: #Otherwise we'll get a division by 0 error
                                        slopeNJ=(ThickWalls[twpid1][4]-ThickWalls[twpid2][4])/(ThickWalls[twpid1][3]-ThickWalls[twpid2][3]) #slope of the line connecting the new junction nodes neighbouring the new wall
                                    else:
                                        slopeNJ=inf
                                    x0=ThickWallNode[3]
                                    y0=ThickWallNode[4]
                                    z0=20.0
                                    #Calculate the horizontal distance between XY0 and the cell center, compare it with the distance between the mean position of the new junctions. If the latter is closer to the cell center, it becomes the new XY0 to make sur the disk is visible
                                    xC=position[int(ThickWallNode[2])][0]
                                    yC=position[int(ThickWallNode[2])][1]
                                    xNJ=(ThickWalls[twpid1][3]+ThickWalls[twpid2][3])/2.0
                                    yNJ=(ThickWalls[twpid1][4]+ThickWalls[twpid2][4])/2.0
                                    if sqrt(square(x0-xC)+square(y0-yC)) > sqrt(square(xNJ-xC)+square(yNJ-yC)):
                                        x0=xNJ
                                        y0=yNJ
                                    for i in range(12):
                                        x=x0+cos(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                        y=y0+sin(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                        z=z0+radiusPlasmodesm_disp*sin(int(i)*pi/6.0)
                                        myfile.write(str(x) + " " + str(y) + " " + str(z) + " \n")
                            else:
                                break #interrupts the for loop in case we reached the new junction nodes
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(len(ThickWalls)+len(PlasmodesmFlowDensity)) + " " + str(len(ThickWalls)*5+len(PlasmodesmFlowDensity)*13) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[1]>=coro: #wall that is a junction
                                myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[5])) + " " + str(int(ThickWallNode[5])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n")
                                myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[6])) + " " + str(int(ThickWallNode[6])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n")
                        for i in range(len(PlasmodesmFlowDensity)):
                            if PlasmodesmFlowDensity[i]==0:
                                myfile.write("12 " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " \n")
                            else:
                                myfile.write("12 " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+1+len(ThickWalls)*2) + " " + str(i*12+2+len(ThickWalls)*2) + " " + str(i*12+3+len(ThickWalls)*2) + " " + str(i*12+4+len(ThickWalls)*2) + " " + str(i*12+5+len(ThickWalls)*2) + " " + str(i*12+6+len(ThickWalls)*2) + " " + str(i*12+7+len(ThickWalls)*2) + " " + str(i*12+8+len(ThickWalls)*2) + " " + str(i*12+9+len(ThickWalls)*2) + " " + str(i*12+10+len(ThickWalls)*2) + " " + str(i*12+11+len(ThickWalls)*2) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(len(ThickWalls)+len(PlasmodesmFlowDensity)) + " \n")
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[1]>=coro: #wall that is a junction
                                myfile.write("9 \n") #Quad cell type
                                myfile.write("9 \n") #Quad cell type
                        for i in range(len(PlasmodesmFlowDensity)):
                            myfile.write("7 \n") #Polygon cell type 
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(ThickWalls)*2+len(PlasmodesmFlowDensity)*12) + " \n")
                        myfile.write("SCALARS TM_n_PD_flux_(m/s) float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[0]<len(MembraneFlowDensity):
                                myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                            else:
                                myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[0]<len(MembraneFlowDensity):
                                myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                            else:
                                myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                        for i in range(len(PlasmodesmFlowDensity)):
                            for j in range(12):
                                myfile.write(str(float(PlasmodesmFlowDensity[i])/sperd/cmperm) + " \n")
                    myfile.close()
                    text_file.close()
        
    
    #write down kr_tot and Uptake distributions in matrices
    iBarriers=-1
    for Barriers in Barrierrange:
        Maturity=int(Barriers.get("value"))
        iBarriers+=1
        text_file = open(newpath+"Macro_prop_"+str(Maturity)+".txt", "w")
        with open(newpath+"Macro_prop_"+str(Maturity)+".txt", "a") as myfile:
            myfile.write("Macroscopic root radial hydraulic properties, apoplastic barrier "+str(Maturity)+" \n")
            myfile.write("\n")
            myfile.write(str(Nscenarios-1)+" scenarios \n")
            myfile.write("\n")
            myfile.write("Cross-section height: "+str(height*1.0E-04)+" cm \n")
            myfile.write("\n")
            myfile.write("Cross-section perimeter: "+str(perimeter[0])+" cm \n")
            myfile.write("\n")
            myfile.write("Cross-section radial conductivity: "+str(kr_tot[iBarriers][0])+" cm/hPa/d \n")
            myfile.write("\n")
            myfile.write("Number of radial discretization boxes: \n")
            r_discret_txt=' '.join(map(str, r_discret.T))
            myfile.write(r_discret_txt[1:21]+" \n")
            myfile.write("\n")
            myfile.write("Radial distance from stele centre (microns): \n")
            for j in Layer_dist2:
                myfile.write(str(float(j))+" \n")
            myfile.write("\n")
            myfile.write("Standard Transmembrane uptake Fractions (%): \n")
            for j in range(int(r_discret[0])):
                myfile.write(str(STFlayer_plus[j][iBarriers]*100)+" \n")
            myfile.write("\n")
            myfile.write("Standard Transmembrane release Fractions (%): \n")
            for j in range(int(r_discret[0])):
                myfile.write(str(STFlayer_minus[j][iBarriers]*100)+" \n")
            for i in range(1,Nscenarios):
                myfile.write("\n")
                myfile.write("\n")
                myfile.write("Scenario "+str(i)+" \n")
                myfile.write("\n")
                myfile.write("h_x: "+str(Psi_xyl[0][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("h_s: "+str(Psi_soil[0][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_x: "+str(Os_xyl[0][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_s: "+str(Os_soil[0][i])+" to "+str(Os_soil[0][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("Xcontact: 999999 microns \n")
                myfile.write("\n")
                myfile.write("kw: "+str(kw)+" cm^2/hPa/d \n")
                myfile.write("\n")
                myfile.write("Kpl: "+str(Kpl)+" cm^3/hPa/d \n")
                myfile.write("\n")
                myfile.write("kAQP: "+str(kaqp_cortex)+" cm/hPa/d \n")
                myfile.write("\n")
                myfile.write("s_hetero: 0 \n")
                myfile.write("\n")
                myfile.write("s_factor: "+str(s_factor[0][count])+" \n")
                myfile.write("\n")
                myfile.write("Os_hetero: 0 \n")
                myfile.write("\n")
                myfile.write("Os_cortex: "+str(Os_cortex[0][count])+" hPa \n")
                myfile.write("\n")
                myfile.write("q_tot: "+str(Q_tot[iBarriers][i]/height/1.0E-04)+" cm^2/d \n")
                myfile.write("\n")
                myfile.write("Stele, cortex, and epidermis uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(UptakeLayer_plus[j][iBarriers][i])+" \n")
                myfile.write("\n")
                myfile.write("Stele, cortex, and epidermis release distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(UptakeLayer_minus[j][iBarriers][i])+" \n")
                myfile.write("\n")
                myfile.write("Xylem uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(abs(Q_xyl_layer[j][iBarriers][i]))+" \n")
                myfile.write("\n")
                myfile.write("Cell layers pressure potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(PsiCellLayer[j][iBarriers][i])+" \n")
                myfile.write("\n")
                myfile.write("Cell layers osmotic potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(OsCellLayer[j][i])+" \n")
                myfile.write("\n")
                myfile.write("Wall layers pressure potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(PsiWallLayer[j][iBarriers][i]/NWallLayer[j][iBarriers][i])+" \n")
                myfile.write("\n")
                myfile.write("Wall layers osmotic potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(OsWallLayer[j][i])+" \n")
        myfile.close()
        text_file.close()
    
    
    #text_file = open(newpath+"Cortical_cell_perimeters.txt", "w")
    #with open(newpath+"Cortical_cell_perimeters.txt", "a") as myfile:
    #    for j in range(10):
    #        myfile.write(str(cortex_cellperimeters[j][:])+" \n")
    #myfile.close()
    #text_file.close()
