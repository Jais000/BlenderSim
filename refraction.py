import bpy
import bmesh
import mathutils
from bpy.app.handlers import persistent
import numpy as np
import math
from numpy import sqrt, sin, cos, pi,power
from cmath import atan
from bpy.types import SpaceView3D
from mathutils.bvhtree import BVHTree
import cmath
import os 
import sys
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir)
from Magnetism import MagModal 
coord = {}
dlocs = {}
#coord[frame] = []
currentTime = 0
projMap ={}
frameselect = 0 
currFrame = 0 
currentFrame = 0
lengths = []
T =.001
pastHits = {}
pastHitsMap = {}
progression = {}
origin = {}
framecount = 0
absdistrec = {}
dirs = {}
dirsMap = {}
flag = {}
lowestDist = {}
lowestDistMap = {}
originDir = {}
loc = {}
originDirMap = {}
locMap = {}
originLoc = {}
originLocMap = {}
project = {}
fsegs = []
source = [] 
absdistrecMap = {}
global finish
present = 0 
global calls 
global count
def createPlane(name,location):
    bpy.ops.mesh.primitive_plane_add(
            calc_uvs=True,
            enter_editmode=False,
            align='WORLD',
            location=location,
            rotation=(0, 0, 0),
            scale=(0, 0, 0))
    bpy.data.objects[bpy.context.selected_objects[0].name].name = name

def transmit(n,i,n_1,n_2):
    n_h = mathutils.Vector((abs(n[0]),abs(n[1]),abs(n[2])))
    
    i_h = mathutils.Vector((abs(i[0]),abs(i[1]),abs(i[2])))
    ni = n.cross(i)
    ni_h = n_h.cross(i_h)
    mu = n_1/n_2
    b = mathutils.Vector(((mu * i) + (mu * ((-n).dot(i)) - cmath.sqrt(1- (mu**2*(1 - ((-n).dot(i))**2)) ) )*n)) 
    ######################print(b)
    if math.isnan(b[0]):
        b = i.reflect(n)    
    return b

name = "RaySourcePlane"


def raycast(loc,dir):
        ###################################print("DIR",dir)
        dg = bpy.context.evaluated_depsgraph_get()
        RESULT, LOCATION, NORMAL, INDEX, OBJECT, MATRIX = bpy.data.scenes[0].ray_cast(dg,loc, dir)
        return [RESULT, LOCATION, NORMAL, INDEX, OBJECT, MATRIX,True]
   




def coords(target,n_1,n_2,S,prevBeamStatus,sourceN,length,frame = -1,delta=True):
     
    #coords(target,1.000293,lengthDict[n],S,dirDict[n],prevBeamStatus[n],str(dirs)+"_"+str(F),length,F,delta)   
    
    global originDirMap
    global pastHits
    global pastHitsMap
    global originLocMap
    global absdistrecMap
    global absdistrec
    global lowestDist
    global lowestDistMap
    global flag
    global T
    global progression
    global calls
    global finish
    global project
    global projMap
    global origin
    global originDir
    global dlocs 
    global coord
    global dirs
    global source
    global originLoc
    global loc
    global locMap 
    global fsegs
    global lengths
    map = {}
    orgmap = {}

    #####################print("LETSSEE",originLoc)

    #######################print("TESTNONANIM",sourceN)
    for n,i in enumerate(source):
        #########################print("NAME",source)
        initLs = [bpy.data.objects[i.name]["Minimal Wavelength"] for i in source] 
        finalLs = [bpy.data.objects[i.name]["Maximal Wavelength"] for i in source]
        lengths = []
        N = [bpy.data.objects[i.name]["Beams per ray"] for i in source ]

       
               
        
        for i in zip(N,zip(initLs,finalLs)):
                        
            lengths.append(np.linspace(i[1][0],i[1][1],i[0]))
            
            

    
    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    for l in lengths[m]:
                        absdistrec[f"{m}_{f}"][l]
        else:
            for m,src in enumerate(source):
                for l in lengths[m]:
                    absdistrec[m][l]      
    except:
        absdistrecMap = {}
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,srcs in enumerate(source):
                ########################print("L",lengths)
                    for l in [k for k in lengths[m]]:
                        absdistrecMap = {}
                    ########################print("l",l)
                        absdistrecMap[l] = 150000
                    absdistrec[f"{m}_{f}"] = absdistrecMap
        else:
            for m,srcs in enumerate(source):
                #############print("KEYERROR",lengths)
                for l in [k for k in lengths[m]]:
                    
                    ########################print("l",l)
                    absdistrecMap[l] = 150000
                    ######################print("L",l)
                absdistrec[m] = absdistrecMap                    
    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    for l in lengths[m]:
                        pastHits[f"{m}_{f}"][l]
        else:
            for m,src in enumerate(source):
                for l in lengths[m]:
                    pastHits[m][l]
    except:
        pastHitsMap = {}
        pastHits = {}
        pastHits[sourceN] = {}
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs:
                for m,src in enumerate(source):            
                    for l in [k for k in lengths[m]]:
                        pastHitsMap[l] = []
                        ##############print("Length",l)
                    pastHits[f'{m}_{f}'] = pastHitsMap
        else:
            for m,src in enumerate(source):            
                for l in [k for k in lengths[m]]:
                    pastHitsMap[l] = []
                    #############print("Length",l)
                pastHits[m] = pastHitsMap
            ##############print(pastHits.keys())
    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    for l in lengths[m]:
                        origin[f"{m}_{f}"][l] 
        else:
            for m,src in enumerate(source):
                for l in lengths[m]:
                    origin[m][l]       
    except:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    orgmap = {}
                    for l in lengths[m]:
                        ###########print(getMatrixData(source)[0])
                        orgmap[l] = getMatrixData(source)[0][m] 
                    origin[f"{m}_{f}"] = orgmap
        else:
            for m,src in enumerate(source):
                orgmap = {}
                #######print(lengths)
                for l in lengths[m]:
                    #######print(m,l)
                    ###########print(getMatrixData(source)[0])
                    orgmap[l] = getMatrixData(source)[0][m] 
                origin[m] = orgmap
    ########print(origin)
    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    dirs[f"{m}_{f}"] 
        else:
            for m,src in enumerate(source):
                dirs[m]       
    except:    
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs:
                for m,src in enumerate(source):
                    dirsMap = {}
                    for l in lengths[m]:
                       
                        dirsMap[l] = getMatrixData(source)[1][m] 
                    dirs[f"{m}_{f}"] = dirsMap
        else: 
            for m,src in enumerate(source):
                dirsMap = {}
                for l in lengths[m]:
                   
                    dirsMap[l] = getMatrixData(source)[1][m] 
                dirs[m] = dirsMap 


    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    originDir[f"{m}_{f}"] 
        else:
            for m,src in enumerate(source):
                originDir[m]      
    except:
        originDirMap = {}
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,srcs in enumerate(source):
                    originDirMap = {}
                    for l in [k for k in lengths[m]]:
                        
                        absdistrecMap[l] = dir[m][l]
                    originDir[f"{m}_{f}"] = originDirMap
        else:
            for m,srcs in enumerate(source):
                originDirMap = {}
                for l in [k for k in lengths[m]]:
                    #print(dirs)
                    originDirMap[l] = dirs[m][l]
                originDir[m] = originDirMap                                              
    
    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    project[f"{m}_{f}"] 
        else:
            for m,src in enumerate(source):
                project[m]   
    except:    
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs:
                for m,src in enumerate(source):
                    projMap = {}
                    for l in lengths[m]:
                        
                        projMap[l] = 0 
                    project[f"{m}_{f}"] = projMap
        else: 
            for m,src in enumerate(source):
                projMap = {}
                for l in lengths[m]:
                    
                    projMap[l] = 0
                project[m] = projMap          

    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    coord[f"{m}_{f}"] 
        else:
            for m,src in enumerate(source):
                coord[m]
    except:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs:
                for m,src in enumerate(source):
                    map = {}
                    for l in lengths[m]:
                        map[l] = [origin[m][l]] 
                    coord[f"{m}_{f}"] = map
        else: 
            for m,src in enumerate(source):
                map = {}
                for l in lengths[m]:
                    print(origin)
                    map[l] = [origin[m][l]] 
                coord[m] = map     
                   

    threshold =140

    meterspercall =  bpy.data.scenes["Scene"].options.speed * T
    newDirs = {}
    newlocs = {}
    convertedDirs ={}

    progMap = {}

    try:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs: 
                for m,src in enumerate(source):
                    progression[f"{m}_{f}"] 
        else:
            for m,src in enumerate(source):
                progression[m]      
    except:
        if bpy.data.scenes["Scene"].options.animating == True:
            for f in fsegs:
                for m,src in enumerate(source):
                    progMap = {} 
                    for l in lengths[m]:
                        progMap[l] = -1 
                    progression[f"{m}_{f}"] = progMap 
        else: 
            for m,src in enumerate(source):
                progMap = {} 
                for l in lengths[m]:
                    progMap[l] = -1 
                progression[m] = progMap 
    #####print(project)
    if project[sourceN][length] == 0: 
        flag = {}
        locMap = {}
        flagMap = {} 
        dirMap ={}
        map = {}
        
        try:
            if bpy.data.scenes["Scene"].options.animating == True:
                for f in fsegs: 
                    for m,src in enumerate(source):
                        flag[f"{m}_{f}"] 
            else:
                for m,src in enumerate(source):
                    flag[m]      
        
        except:
            if bpy.data.scenes["Scene"].options.animating == True:
                for f in fsegs:
                    for m,src in enumerate(source):
                        flagMap = {}
                        for l in lengths[m]:
                            flagMap[l] = 0
                        flag[f"{m}_{f}"] = flagMap 
            else: 
                for m,src in enumerate(source):
                    flagMap = {} 
                    for l in lengths[m]:
                        flagMap[l] = 0
                    flag[m] = flagMap                
      
        

        try: 
            if bpy.data.scenes["Scene"].options.animating == True:
                for f in fsegs: 
                    for m,src in enumerate(source):
                        coord[f"{m}_{f}"] 
            else:
                for m,src in enumerate(source):
                    coord[m]      
        except:
            if progression[m][l]<0:
                if bpy.data.scenes["Scene"].options.animating == True:
                    for f in fsegs:
                        for m,src in enumerate(source):
                            map = {}
                            for l in lengths[m]:
                                map[l] = [origin[m][l]]
                            coord[f"{m}_{f}"] = map 
            else: 
                if progression[m][l]<0:
                    for m,src in enumerate(source):
                        map = {} 
                        for l in lengths[m]:
                            map[l] = [origin[m][l]]
                        coord[m] = map                   

        
        if progression[sourceN][length]<0:
            lowestDistMap[length] = threshold
            lowestDist[sourceN] = lowestDistMap
            loc = {}
            coord = {}
            
            if bpy.data.scenes["Scene"].options.animating == True:
                for f in fsegs:
                    for m,srcs in enumerate(sources): 
                        for l in [k for k in lengths[m]]:
                            locMap[l] = origin[m][l]
                        loc[f"{m}_{f}"] = locMap           
            else:
                for m,srcs in enumerate(source):  
                    for l in [k for k in lengths[m]]:
                        print(origin)
                        locMap[l] = origin[m][l]
                    loc[m] = locMap     

            newDirs = {}
            convertedDirs = {}
            newlocs = {}
            try:
                if bpy.data.scenes["Scene"].options.animating == True:
                    for f in fsegs: 
                        for m,src in enumerate(source):
                            coord[f"{m}_{f}"] 
                else:
                    for m,src in enumerate(source):
                        coord[m]
            except:
                if bpy.data.scenes["Scene"].options.animating == True:
                    for f in fsegs:
                        for m,src in enumerate(source):
                            map = {}
                            for l in lengths[m]:
                                map[l] = [origin[m][l]] 
                            coord[f"{m}_{f}"] = map
                else: 
                    for m,src in enumerate(source):
                        map = {}
                        for l in lengths[m]:
                            map[l] = [origin[m][l]] 
                        coord[m] = map     
        #####################print("ARRAY",origin)
            try:
            
                if bpy.data.scenes["Scene"].options.animating == True:
                    for f in fsegs: 
                        for m,src in enumerate(source):
                            originLoc[f"{m}_{f}"] 
                else:
                    for m,src in enumerate(source):
                        originLoc[m]      
            except:
                originLocMap = {}
                originLoc = {}
                if bpy.data.scenes["Scene"].options.animating == True:            
                    for f in fsegs: 
                        for m,srcs in enumerate(source):
                            for l in [k for k in lengths[m]]:
                            
                                originLocMap[l] = loc[m][l]
                            originLoc[str(m)+"_"+str(f)] = originLocMap
                else:
                    for m,srcs in enumerate(source):
                        #############print("LENGTHS",lengths)
                        for l in [k for k in lengths[m]]:
                            ############print(loc)
                            originLocMap[l] = loc[m][l]
                        originLoc[m] = originLocMap          
        ###print(coord)
        for k,length in enumerate(dirs[sourceN].keys()):
            ray = []
            try:

                ray=raycast(loc[sourceN][length],dirs[sourceN][length])
            except:

                ray=raycast(loc[sourceN][length][0],dirs[sourceN][length])

            if ray[0] == True:
                if pastHits[sourceN][length] == []:
                    pastHits[sourceN][length] = [ray]
                dist = loc[sourceN][length] - pastHits[sourceN][length][project[sourceN][length]][1]
                absdist = sqrt(dist[0]**2 + dist[1]**2 + dist[2]**2)
                absdistrec[sourceN][length] = min(absdistrec[sourceN][length],absdist)
                if absdist > absdistrec[sourceN][length] or bpy.data.scenes["Scene"].options.animating == False:
                    
                    if ray[4].name in [i.name for i in target]:
                        normal = ray[2]
                        newloc = ray[1]
                        dir_ = dirs[sourceN][length]
                        ior = n_2[length][ray[4].name]
                        pastHits[sourceN][length].append(ray)
                        loc[sourceN][length] = ray[1]
                        dirs[sourceN][length] = transmit(normal,dir_,n_1,ior)
                        origin[sourceN][length] = ray[1]
                        if bpy.data.scenes["Scene"].options.animating==True:
                            flag[sourceN][length] = True
                            loc[sourceN][length] = newloc
                            
                            coord[sourceN][length].append(newloc)  
                            
                            progression[sourceN][length] = 0
                        else: 
                            #######################print("SOURCE",sourceN)
                            flag[sourceN][length] = True
                            loc[sourceN][length] = newloc
                            #####print("COORDBROKE",coord)
                            coord[sourceN][length].append(newloc)
                            progression[sourceN][length] = 0
                            
                        newlocs[length] = newloc
                        #############################print("IN",newlocs[length])
                    else:
                        if bpy.data.scenes["Scene"].options.animating == True:
                            newloc = ray[1]
                            newDirs[length] = dirs[sourceN][length]
                            newlocs[length] = newloc
                            coord[sourceN][length].append(ray[1])
                            flag[sourceN][length] = True
                        else:
                            coord[sourceN][length].append(ray[1])
                            newloc = ray[1]
                            newDirs[length] = dirs[sourceN][length]
                            newlocs[length] = newloc
                else:
                    counting = 0 
                    for i in pastHits[sourceN][length]:
                        counting += 1                        
                    if not pastHits[sourceN][length][counting-1]  == ray:
                        pastHits[sourceN][length].append(ray)
                    fram = bpy.context.scene.frame_current
                    dloc = loc[sourceN][length]+ dirs[sourceN][length]* bpy.data.scenes["Scene"].options.speed*T
                    loc[sourceN][length] = dloc
                    progression[sourceN][length] += 1
                    newDirs[length] = dirs[sourceN][length]
                    ########################print("flag",fram, progression[sourceN][length],coord[sourceN][length])
                    flag[sourceN][length] = False
            else:
                dist = loc[sourceN][length]-ray[1]
                absdist = sqrt(dist[0]**2 + dist[1]**2 + dist[2]**2)
                ########################print("FAIL",absdist)
                if bpy.data.scenes["Scene"].options.animating==True:
                    # originDir[sourceN.split("_")[0]][length]*S
                    
                    coord[sourceN][length].append(coord[sourceN][length][0] + originDir[int(sourceN.split("_")[0])][length]*S) 
                else:
                    
                    
                    try:
                        ###################print("COORD",coord)
                        coord[sourceN][length].append(coord[sourceN][length][0] + originDir[length]*S)
                        #######################print("COORD ERROR",coord[sourceN][length][project[sourceN][length]+1])
                        newloc = coord[sourceN][length][project+1] + originDir[length]*S
                    except:
                        ###################print("COORD ERROR",coord[sourceN][length][0] ,originDir[sourceN][length]*S)
                        coord[sourceN][length].append(coord[sourceN][length][0] + originDir[sourceN][length]*S)
                        
                        newloc = coord[sourceN][length][project[sourceN][length]+1] + originDir[sourceN][length]*S                            
                newDirs[length]=dirs[sourceN][length]
                newlocs[length]=loc[sourceN][length]
                    
        ################################print("TRUTH",flag)
        #casee = flag["0_0"][length]
       
        if flag[sourceN][length] == True:
            #######################print("NEW COORD",sourceN)
            project[sourceN][length] += 1 
            progression[sourceN][length] = -1
            deletor = []
            flag[sourceN][length] == False
            coords(target,n_1,n_2,S,prevBeamStatus,sourceN ,length,frame )
    
    elif project[sourceN][length] <= bpy.data.scenes["Scene"].options.RecursiveDepth and project[sourceN][length] > 0:
        #######################print("project",sourceN,project)
        #########################print("COOORDS",coord)

        newDirs = {}
        newlocs = {}
        locmap = {}
        first = ""
        second = ""
        mx = 0
        ################################print("COOORDDS",coord)

            #################################print("MAXXES", first,second)
        if bpy.data.scenes["Scene"].options.animating == True  or bpy.data.scenes["Scene"].options.animating == False:
            counting = 0 
            for i in coord[sourceN][length]:
                counting += 1
                ########################print(counting)
            #################################print("COORD first second ",coord[first][second],len(coord[first][second]))
            if project == counting:
                progMap = {}
                for n,wavelength in enumerate(dir):
                    progMap[wavelength] = 0
                    #locmap[wavelength] = [] 
                    #################################print("HASHI",coord[[i for i in coord][sourceN]][wavelength])
                    #coord[str([i for i in coord][sourceN].split("_")[0])+"_"+str(frame)][wavelength] = locmap
                progression[sourceN] = progMap
            try:
                 dlocs[sourceN][length]
            except:
                locmap[length] = [] 
                dlocs[sourceN] = locmap
            for k,length in enumerate([length]):
                ray = raycast(loc[sourceN][length]+.0001*dirs[sourceN][length],dirs[sourceN][length])
                #######################print("PROG",progression[sourceN][length])

                if ray[0] == True:
                    ################print(pastHits)
                    if pastHits[sourceN][length] == []:
                        pastHits[sourceN][length] = [ray]
                    
                    #######################print("PAST",project[sourceN][length],"PAST HIT:",pastHits[sourceN][length][1])
                    ##############print("LENGTH",pastHits[sourceN][length][project[sourceN][length]-1][1])
                    dist = loc[sourceN][length] - pastHits[sourceN][length][project[sourceN][length]-1][1]
                    absdist = sqrt(dist[0]**2 + dist[1]**2 + dist[2]**2)
                    absdistrec[sourceN][length] = min(absdistrec[sourceN][length],absdist)
                    if absdist > absdistrec[sourceN][length] or bpy.data.scenes["Scene"].options.animating == False:
                        
                        if ray[4].name in [i.name for i in target]:
                            ##################print("COOOORDS",project[sourceN][length],coord)
                            pastHits[sourceN][length].append(ray)
                            normal = ray[2]
                            dir_ = dirs[sourceN][length]
                            ior = n_2[length][ray[4].name]
                            #bpy.ops.object.empty_add(location=(dirs[sourceN][length]]))
                            if dir_.dot(normal)>0:
                               newDir = transmit(-normal,dir_,ior,n_1)
                            else:
                                newDir = transmit(normal,dir_,n_1,ior)
                            dirs[sourceN][length] = newDir
                            loc[sourceN][length] = ray[1]
                            if bpy.data.scenes["Scene"].options.animating == True:
                                
                                ########################print("PASS",coord)
                                coord[sourceN][length].append(ray[1])
                                newDirs[length] = newDir
                                newloc = ray[1]
                                newlocs[length] = newloc
                                project[sourceN][length] += 1
                                flag[sourceN][length] == True
                            else:
                                coord[sourceN][length].append(ray[1])
                                newDirs[length] = newDir
                                newloc = ray[1]
                                newlocs[length] = newloc
                                project[sourceN][length] += 1
                                flag[sourceN][length] == True                                

                        else:
                            if bpy.data.scenes["Scene"].options.animating == True:
                                coord[sourceN][length].append(ray[1])
                            #coord[str(sourceN)+"_"+str(frame)][length].append(ray[1])
                            else:
                                coord[sourceN][length].append(ray[1])
                            #coord[sourceN][length].append(ray[1])
                            newloc = ray[1]
                            newDirs[length] = dirs[sourceN][length]
                            newlocs[length] = newloc
                      
                    else:
                        counting = 0 
                        for i in pastHits[sourceN][length]:
                            counting += 1                        
                        if not pastHits[sourceN][length][counting-1]  == ray:
                            pastHits[sourceN][length].append(ray)
                        dloc = loc[sourceN][length] + dirs[sourceN][length]* bpy.data.scenes["Scene"].options.speed*T
                        loc[sourceN][length] = dloc
                        progression[sourceN][length] += 1
                        newDirs[length] = dirs[sourceN][length]
                        #########################print(flag)
                        flag[sourceN][length] = False
                        
                else:
                    newDirs[length] = dirs[sourceN][length]
                    newlocs[length] = loc[sourceN][length]+S*dirs[sourceN][length]
                    if prevBeamStatus[length] == True:
                        if bpy.data.scenes["Scene"].options.animating == True:
                            counting = 0
                            for i in coord[sourceN][length]:
                                counting += 1
                            where = counting-1
                            coord[sourceN][length].append(coord[sourceN][length][where] +S*dirs[sourceN][length]  )
                            #coord[str(sourceN)+"_"+str(frame)][length].append(coord[str(sourceN)+"_"+str(frame)][length][len(coord[str(sourceN)+"_"+str(frame)][length])-1] + S*dir[length])
                        #bpy.ops.object.empty_add(location=(dir[length]))
                        else:
                            counting = 0
                            for i in coord[sourceN][length]:
                                counting += 1
                            coord[sourceN][length].append(coord[sourceN][length][counting-1] +S*dirs[sourceN][length]  )
                            #coord[sourceN][length].append(coord[sourceN][length][len(coord[sourceN][length])-1] + S*dir[length])
                        ####################################print("COOORD",coord)
                        prevBeamStatus[length] = False
        
        if flag[sourceN][length] ==True :
            progression[sourceN][length] == -1
            project[sourceN][length] += 1
            flag[sourceN][length] == False
            coords(target,n_1,n_2,S,prevBeamStatus,sourceN,length,frame)
    else:
        #coord = {}
        ##################print("FINAL",coord)
        progression = {}
        project = {}
        dirs = {}
        flag ={}
    
def array(num,obj,x=-1):
    mod = bpy.data.objects[obj.name].modifiers.new(name='array', type='ARRAY')
    mod.use_constant_offset = True
    mod.use_relative_offset = False
    mod.constant_offset_displace[0] = .1
    mod.fit_type = 'FIT_LENGTH'
    mod.fit_length = bpy.data.scenes["Scene"].options.photonLength
    deform = bpy.ops.object.modifier_add(type='CURVE')
   

    #bpy.context.collection.objects.link(obj)
    #################################print("NUMBER",num,x)
    if bpy.data.scenes["Scene"].options.animating == True:
        index = obj.name.split(".")
        
        if len(index)>1:
            #################################print("INDEX")
            #mod.curve = bpy.data.objects["LightPath"+str(num)+"_"+str(x)+"."+index[1]]
            #deform = bpy.ops.object.modifier_add(type='CURVE')
          
            bpy.context.object.modifiers["Curve"].object = bpy.data.objects["LightPath"+str(num)+"_"+str(x)+"."+index[1]]
        else:
            ###########################print("COOORD",coord)
            obj.modifiers["Curve"].object = bpy.data.objects["LightPath"+str(num)+"_"+str(x)]
            #################################print("LightPath"+str(num)+"_"+str(x))
        
    else:    
        mod.curve = bpy.data.objects["LightPath"+str(num)]
        deform = bpy.ops.object.modifier_add(type='CURVE')
        bpy.context.object.modifiers["Curve"].object = bpy.data.objects["LightPath"+str(num)]
     
    
def createRenderObj(N,i,data,number=-1):
    radius = bpy.data.scenes["Scene"].options.photonRadius
    bpy.ops.mesh.primitive_cylinder_add(radius=radius, depth=.1, enter_editmode=True, align='WORLD', rotation=(pi/2,0,0),location=(0, 0, 0), scale=(1, 1, 1))
    bpy.data.objects[bpy.context.selected_objects[0].name].name = "render-"+str(i)+"-"+str(number)
    bpy.ops.transform.rotate(value=1.5708, orient_axis='Z', orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', constraint_axis=(False, False, True), mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    bpy.ops.object.editmode_toggle()
    lightRender = bpy.context.selected_objects[0]
    sfx = lightRender.name.split(".")
    if len(sfx)>1:
        sfx = "."+ sfx[1]
    else:
        sfx = ""
    if bpy.data.scenes["Scene"].options.animating==True:
        array(i,lightRender,number)
        finalFrame = bpy.data.scenes["Scene"].options.endTime
        animate(lightRender,0,4000,int(number),finalFrame)
    

    else:
        array(i,lightRender)
        animate(lightRender,0,1000)


def animate(obj,i,f,s=0,e=1):

  
    bpy.ops.object.editmode_toggle()
    bpy.context.scene.frame_current = int(s)
    bpy.ops.transform.translate(value=(i, -0, -0), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', constraint_axis=(True, False, False), mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    bpy.ops.anim.insert_keyframe_animall()
    fps = bpy.context.scene.render.fps
    speed = bpy.data.scenes["Scene"].options.speed
    bpy.context.scene.frame_current = int(e)
    bpy.ops.transform.translate(value=(speed*((e-s)/fps), 0, 0), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', constraint_axis=(True, False, False), mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    bpy.ops.anim.insert_keyframe_animall()
    bpy.ops.object.editmode_toggle()
    bpy.context.scene.frame_current = int(s)
    
def dispersion(dir,target,loc,S,an,nd,N,sources,orig,iors,F=-1,possible = [-1],delta=True  ):
    
    poss = []
    initLs = [bpy.data.objects[i.name]["Minimal Wavelength"] for i in sources] 
    finalLs = [bpy.data.objects[i.name]["Maximal Wavelength"] for i in sources]
    lengths = []
    global origin
    global originDir
    prevBeamStatus = {}
    for i in zip(N,zip(initLs,finalLs)):
        lengths.append(np.linspace(i[1][0],i[1][1],i[0]))
    lengthDict={}
    dict = {}
    progression = {}
    dirDict = {}
    locDict = {}
    dxs = []
    B = []
    A = []
    
    for n,i in enumerate(an):
        ###################################print("IORS",iors,orig)
        B.append(((iors[n]-1)/i)*0.52)
        A.append(iors[n] - (B[n]*2.897))
    for v,j in enumerate(lengths):#direction
        dict = {}
        for o,k in enumerate(j):#lengths
            map = {}
            for n,prism in enumerate(target):
                map[prism.name]=A[n]+(B[n]/(k**2))
            dict[k] = map
        lengthDict[v] = dict
    
    for n,i in enumerate(sources):
        
        #dirDict[n]={j:dir[n] for j in lengths[n]}
        prevBeamStatus[n] = {j:True for j in lengths[n]}
        #########################print("lengths",prevBeamStatus)
        ###################################print("orig",orig)
    if bpy.data.scenes["Scene"].options.animating == True:
        if orig == True:
                for n,dirss in enumerate(lengthDict):
                    ################################print('dirs',dirs)
                    for length in lengthDict[dirss]:
                        try:
                        ###################################print("LENGTHS",lengths,[t for t in [h.values() for h in [u for u in [s for s in loc.values()]]][0]])
                            #########################print("STARTINGH",lengthDict,prevBeamStatus)
                            origin[dirss]={j:[t for t in [h.values() for h in [u for u in [s for s in loc.values()]]][0]][n] for j in lengths[n]}
                            
                            originDir[dirss] = {j:dir[n] for j in lengths[n]}
                        except:
                        ###################################print("error",loc,lengths)
                            ##############################print("BROKE",origin, loc,lengths)
                            try:
                                #########################print("ERROR",dir)
                                origin[dirss]={j:loc[dirss] for j in lengths[dirss]}
                                originDir[dirss] = {j:dir[dirss] for j in lengths[dirss]}
                            except:
                                #########################print("ERROR",loc,dir)
                                origin[dirss]={j:loc[n] for j in lengths[n]}
                                originDir[dirss] = {j:dir[n] for j in lengths[n]}
                        #########################print(lengthDict, dirDict)
                        coords(target,1.000293,lengthDict[n],S,prevBeamStatus[n],str(dirss)+"_"+str(F),length,F,delta)   
        else:
                for n,dirss in enumerate(lengthDict):
                    for length in lengthDict[dirss]:
                        ##############################print("DISP",lengthDict,locDict,dirDict)
                        coords(target,1.000293,lengthDict[n],S,prevBeamStatus[n],str(dirss)+"_"+str(F),length,F,delta)   
    else:
        for n,dirss in enumerate(lengthDict):
                    ################################print('dirs',dirs)
            for length in lengthDict[dirss]:
                try:
                                #########################print("ERROR",dir)
                    origin[dirss]={j:loc[dirss] for j in lengths[dirss]}
                    originDir[dirss] = { j:dir[dirss] for j in lengths[dirss]}
                    
        
                except:
                                #########################print("ERROR",loc,dir)
                    origin[dirss]={j:loc[n] for j in lengths[n]}
                    originDir[dirss] = {j:dir[n] for j in lengths[n]}
                    
                for n,dirs in enumerate(lengthDict):
                    for length in lengthDict[dirs]:
                        locDict[n] = {j:loc[n] for j in lengths[n]}
                        coords(target,1.000293,lengthDict[n],S,prevBeamStatus[n],dirs,length,F,delta)   
    return lengthDict
def render(N,i,number,data,x):
    ###################################print(number)
    if not number == -1:
        object = bpy.data.objects['render'+"-"+str(i)+"-"+str(number)+x]
    else:
        object = bpy.data.objects['render-'+str(i)+x]
    if not number == -1:
        object = bpy.data.objects['render'+"-"+str(i)+"-"+str(number)+x]
    else:
        object = bpy.data.objects['render-'+str(i)+x]
    bpy.context.view_layer.objects.active = object
    object.visible_camera = False
    object.visible_glossy = False
    name = "Light-"+str(data) 
    try:
        material = bpy.data.materials[name]
        bpy.context.view_layer.objects.active.active_material= material
    except:
        material = bpy.data.materials.new(name= name)
        bpy.context.view_layer.objects.active.active_material = material
        material.use_nodes=True
    area = bpy.context.area
    old_type = area.ui_type
    area.ui_type = 'ShaderNodeTree'

    try:
        material.node_tree.nodes.remove(material.node_tree.nodes["Principled BSDF"])
        nodes = material.node_tree.nodes
        for node in nodes:
            if node.type != 'OUTPUT_MATERIAL': # skip the material output node as we'll need it later
                nodes.remove(node) 
    except:
        nodes = material.node_tree.nodes
        for node in nodes:
            if node.type != 'OUTPUT_MATERIAL': # skip the material output node as we'll need it later
                nodes.remove(node) 
        ###################################print("not")
    material_output = material.node_tree.nodes.get('Material Output')
    volume = material.node_tree.nodes.new("ShaderNodeVolumePrincipled")
        #emission.inputs['Strength'].default_value = 5.0
    volume.inputs['Density'].default_value = 0
    volume.inputs['Emission Strength'].default_value = 5  
    ramp = material.node_tree.nodes.new("ShaderNodeValToRGB")
    elements = material.node_tree.nodes["Color Ramp"].color_ramp.elements
    elements[0].color = (0,0,0,1)
    elements[0].position = 0.0
    elements.new(position = 0.745)
    elements.new(position = 0.7)
    elements.new(position = 0.6)
    elements.new(position = 0.5)
    elements.new(position = 0.4)
    elements.new(position = 0.3)
    elements[6].color = (0,0,0,1)
    elements[5].color = (1,0,0.000795,1)
    elements[4].color = (0,.875,0.000083,1)
    elements[3].color = (0,0.152711,0.75,1)
    elements[2].color = (0.5,0.152711,0.363423,1)
    elements[1].color = (0,0,0,1)
    material.node_tree.links.new(material_output.inputs[1], volume.outputs[0])
    material.node_tree.links.new(volume.inputs[0], ramp.outputs[0])
    value = material.node_tree.nodes.new("ShaderNodeValue")
    value.outputs[0].default_value = data
    ###################################print("DAAAAAAATA",data)
    material.node_tree.links.new(ramp.inputs[0],value.outputs[0])
    material.node_tree.links.new(ramp.outputs[0],volume.inputs[7])
    area.type = old_type
    area.ui_type = old_type
bpy.data.scenes['Scene']['prismDict']
bpy.data.scenes['Scene']['AbbesDict']



                            

#N= bpy.data.scenes["Scene"].my_tool.BandDefinition
class Modal(bpy.types.Operator):
    bl_idname = "timed.ray"
    bl_label = "timedRay"
    
    def modal(self, context, event):
        global currFrame
        global animCount
        global fsegs
        global T 
        global coord 
        if event.type in {'RIGHTMOUSE', 'ESC'}:
            self.cancel(context)
            return {'CANCELLED'}

        if event.type == 'TIMER':
            
            sources = [i for i in bpy.data.objects if "Ray Source" in i.values()]
            locs = [i.location for i in sources]
            eulers = [mathutils.Euler((i.rotation_euler[0], i.rotation_euler[1],i.rotation_euler[2]), 'XYZ') for i in sources]
            dirs = [mathutils.Vector((0,0,1)) for i in sources]
            [i.rotate(j) for k,i in enumerate(dirs) for n,j in enumerate(eulers) if k == n]
            loc=mathutils.Vector((-20,-80,100))


            
            prismsN = bpy.data.scenes['Scene']["PrismCount"]

            targets = [bpy.data.scenes["Scene"][f"prism_{i}"] for i in range(prismsN) if bpy.data.scenes["Scene"][f"prism_{i}"] and bpy.data.scenes["Scene"][f"prism_{i}"] in [l for l in bpy.data.scenes["Scene"].objects]]
            IORs = [bpy.data.scenes["Scene"][f"IOR_{i}"] for i in range(len(targets)) if bpy.data.scenes["Scene"][f"prism_{i}"]]
            aNs = [bpy.data.scenes["Scene"][f"aN_{i}"] for i in range(len(targets)) if bpy.data.scenes["Scene"][f"prism_{i}"]]
            defs = [bpy.data.objects[i.name]["Beams per ray"] for i in sources]
            S= 50

            if bpy.data.scenes["Scene"].options.animating == False:
                for i in bpy.context.scene.objects:
                    try:
                        if "Light" in i.name and not len(i.name.split("_"))>1 :
                            objs = bpy.data.objects
                            objs.remove(i, do_unlink=True)
                            bpy.data.curves.remove(i.data)
                        for i in bpy.data.objects:
                            if not i in bpy.context.scene.objects:
                                objs.remove(i, do_unlink=True)
                                bpy.data.curves.remove(i.data)
                         
                    except:
                        pass
                for i in bpy.data.curves:
                    try:
                        if "crv" in i.name and not len(i.name.split("_"))>1 :
                            bpy.data.curves.remove(i) 
                    except:
                        pass
            elif bpy.context.scene.frame_current <= bpy.data.scenes["Scene"].options.endTime:
            
                ########################print("SHOP",coord)
                directions = [i for i in coord]
                wavelengths = [[j for j in coord[i].keys()] for i in directions]
                for n,j in enumerate(wavelengths):
                    for i in bpy.context.scene.objects:
                        try:
                            if "Light" in i.name:
                                if str(n) in i.name: 
                                    objs = bpy.data.objects
                                    objs.remove(i, do_unlink=True)
                        except:
                            if "Light" in i.name:
                                if str(n) in i.name: 
                                    objs = bpy.data.objects
                                    objs.remove(i, do_unlink=True)
                                    bpy.data.curves.remove(i.data)
                for i in bpy.data.curves:
                    if "crv" in i.name:
                        bpy.data.curves.remove(i)
        
                                           
                        
                                    
            if bpy.data.scenes["Scene"].options.animating == True: 
                lengthDict = getDisp(dirs,targets,dlocs,S,aNs,1.7387,defs,sources,IORs,True)
            else:
                lengthDict = getDisp(dirs,targets,locs,S,aNs,1.7387,defs,sources,IORs,False)
        
        
                if bpy.data.scenes["Scene"].options.animating==True:
                    if bpy.context.scene.frame_current<=bpy.data.scenes["Scene"].options.endTime:
                       
                        for k,dir in enumerate([i for i in coord if len(str(i).split("_"))>1]):
                        #if str(bpy.context.scene.frame_current) == str(dir.split("_")[1]): 
                        ###################################print("VALUES",len([k for k in coord[dir].values()][0]))
                                for i,beam in enumerate([coord[dir][i] for i in coord[dir] if len(coord[dir][i])>1]):
                                    ##########################print("beam",beam)
                                    crv = bpy.data.curves.new('crv_'+str(dir), 'CURVE')
                                    crv.dimensions = '3D'
                                    spline = crv.splines.new(type='POLY')
                                    spline.points.add(len(beam)-1) 
                                    try:
                                        for p, new_co in zip(spline.points,[beam]):
                                            ###########################print("new_co",new_co)
                                            p.co = (np.array((new_co.x,new_co.y,new_co.z,1.0)))
                             
                                        obj = bpy.data.objects.new('LightPath'+str(dir.split("_")[0])+"_"+dir.split("_")[1], crv)
                        
                                        bpy.context.scene.collection.objects.link(obj)
                                    except:
                                        for p, new_co in zip(spline.points,beam):
                                            ###########################print("new_co",new_co)
                                            p.co = (np.array((new_co.x,new_co.y,new_co.z,1.0)))
                             
                                        obj = bpy.data.objects.new('LightPath'+str(dir.split("_")[0])+"_"+dir.split("_")[1], crv)
                        
                                        bpy.context.scene.collection.objects.link(obj)
                                        
                            
                           
                
                        if bpy.context.scene.frame_current == bpy.data.scenes["Scene"].options.endTime:
                            if animCount == 0:
                                for i in coord:
                                    if not len(str(i)) == 1 and len(str(i).split("_"))>1:
                                        #################################print("FRAME",i,coord)     
                                        frame = str(i).split("_")[1]
                                    else:
                                        #################################print("FRAME",i)
                                        frame = "0"
                                    for o in coord[i]:
                                    
                                        #################################print("TEST",str(i).split("_")[0],frame)

                                        createRenderObj("N",str(i).split("_")[0],"data",frame)

                            animCount += 1
                            return {'FINISHED'}
                ###################################print(bpy.context.scene.frame_current)
                else:
                    for k,dir in enumerate(coord):        
                        for i,beam in enumerate(coord[dir].values()):
                            crv = bpy.data.curves.new('crv', 'CURVE')
                            crv.dimensions = '3D'
                            spline = crv.splines.new(type='POLY')
                            spline.points.add(len(beam)-1) 
                            for p, new_co in zip(spline.points,  beam):
                                ####################print([i for i in beam])
                                p.co = (np.array((new_co.x,new_co.y,new_co.z,1.0)))
                         
                            obj = bpy.data.objects.new('LightPath'+str(k), crv)
                        
                            bpy.context.scene.collection.objects.link(obj)
                        
                        #createRenderObj("N",k,"data",x)
                currFrame = currFrame + (T * bpy.context.scene.render.fps)
                ##########################print("CUURRR",currFrame)
                bpy.context.scene.frame_current = round(currFrame)
        return {'PASS_THROUGH'}
    
    def execute(self, context):
        global dlocs 
        global currFrame
        currFrame = 0
        dlocs = {}
        global animCount
        animCount = 0 
        global present
        present = 0 
        global currentTime
        currentTime = 0
        global T
        global coord
        global framecount
        coord = {}
        global project
        global progression
        global loc
        global originDir
        global dirs
        project = {}
        loc = {}
        dirs = {}
        coord = {}
        originDir = {}
        progression = {}
        sceneObjs = bpy.context.scene.objects
        allObjs = bpy.data.objects
        framecount = 0
        for l in allObjs:
            try:
                if "render" in l.name:
                    mesh = l.data
                    allObjs.remove(l,do_unlink=True)
                    bpy.data.meshes.remove(mesh) 
                if "LightPath" in l.name:
                    curve =  l.data
                    allObjs.remove(l,do_unlink=True)
                    bpy.data.curves.remove(curve)       
            except:
                pass
        if bpy.data.scenes["Scene"].options.animating == True:
            currentFrame=bpy.data.scenes["Scene"].options.startTime 
            bpy.context.scene.frame_current = currentFrame
        wm = context.window_manager
        self._timer = wm.event_timer_add(T, window=context.window)
        wm.modal_handler_add(self)
        
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
    

def getMatrixData(sources):
    locs = [mathutils.Vector((i.matrix_world[0][3],i.matrix_world[1][3],i.matrix_world[2][3])) for i in sources]
    dirs = [mathutils.Vector((i.matrix_world[0][2],i.matrix_world[1][2],i.matrix_world[2][2])) for i in sources]
    [i.normalize() for i in dirs]
    ######################print("MATRIX DATA",locs)
    return [locs, dirs]
    

def getDisp(dirs,targets,local_locs,S,aNs,nd,defs,sources,iors,delta):
    global source
    sources = [i for i in bpy.data.objects if "Ray Source" in i.values()]
    source = sources
    
    locs = [i.location for i in sources]
    matrixData = getMatrixData(sources)
    eulers = [mathutils.Euler((i.rotation_euler[0], i.rotation_euler[1],i.rotation_euler[2]), 'XYZ') for i in sources]
    dirs = [mathutils.Vector((0,0,1)) for i in sources]
    
    [i.rotate(j) for k,i in enumerate(dirs) for n,j in enumerate(eulers) if k == n]
    loc=mathutils.Vector((-20,-80,100))
  
    prismsN = bpy.data.scenes['Scene']["PrismCount"]

    targets = [bpy.data.scenes["Scene"][f"prism_{i}"] for i in range(prismsN) if bpy.data.scenes["Scene"][f"prism_{i}"]and bpy.data.scenes["Scene"][f"prism_{i}"] in [l for l in bpy.data.scenes["Scene"].objects]]
    if not bpy.context.object:
        current = bpy.data.objects[0].name

    else:
        current = bpy.context.object.name
    for i in targets:
        i.select_get()
        if bpy.context.object in targets:
            bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
            bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='BOUNDS')
    bpy.data.objects[current].select_get()
    global T
    global frameselect
    global fsegs
    global present
    global currentFrame
    global originDir
    global framecount
    global originLoc
    global currentTime
    global project
    dirs = {} 
    global dirsMap
    if bpy.data.scenes["Scene"].options.animating == True: 
        
        fps = bpy.context.scene.render.fps
        pF = bpy.data.scenes["Scene"].options.photonFreq
        currentTime += T 
        #currentFrame += round(currentTime*fps)
        
        
        sF = bpy.data.scenes["Scene"].options.startTime
        fF = bpy.data.scenes["Scene"].options.endTime
        
        segs = np.arange(sF/fps,fF/fps,1/pF)
        fsegs = [round(i*fps) for i in segs]     
        if bpy.context.scene.frame_current in fsegs:
            frameselect = bpy.context.scene.frame_current
        else:
            for o,f in enumerate(fsegs):
                if bpy.context.scene.frame_current > f:
                    try:
                        if bpy.context.scene.frame_current < fsegs[o+1]:
                            bpy.context.scene.frame_current == f
                    except:
                        framseselect = f
                    
            #################################print("FSEGS",fsegs,fsegs[framecount-1])
        if bpy.context.scene.frame_current in fsegs or not dlocs == {}:
            frameselect = bpy.context.scene.frame_current
            ##########################print("CURREt",frameselect)
            t=0
            for o,f in enumerate(fsegs):
                if bpy.context.scene.frame_current > f:
                    try:
                        if frameselect < fsegs[o+1]:
                            t == fsegs[o+1]
                            ##########################print("FRAMESELECT",t)
                    except:
                        ##########################print("EXCEPT",fsegs[o+1],f)
                        t = f
                if bpy.context.scene.frame_current == f:
                    ##########################print("FRAMESELECT",t)
                    t == f
                if bpy.context.scene.frame_current< f:
                    try:
                        if frameselect > fsegs[o-1]:
                            t = fsegs[o-1]
                            ##########################print("FRAMESELECT",fsegs[o+1],f)
                    except:
                    
                        t = fsegs[o+1]
            frameselect = t
            ##########################print("FRAMESELECT",t,fsegs)
            test = False
            ###################################print("RUNNING",dlocs)
            if bpy.data.scenes["Scene"].options.animating == True:
                if not dlocs.values() ==  [{}] and not dlocs == {}:
                    ###################################print("VALS",[i for i in dlocs.values()])
                    for s in progression: 
                        for wavelength in progression[s]:
                            ###################################print("RUNNING2",progression)
                            if not progression[s][wavelength] == -1:
                                test = True
                        if test == True:
                            ###########################print("FRAME",bpy.context.scene.frame_current ,frameselect)
                            ##########################print("GOING",frameselect)
                            return dispersion(matrixData[1],targets,dlocs,S,aNs,nd,defs,sources,False,iors,frameselect,delta)
                        else:
                            #############################print("GETDISP",dlocs)
                            return dispersion(matrixData[1],targets,dlocs,S,aNs,nd,defs,sources,True,iors,frameselect,delta)
                else:
                   
            
                    for n,i in enumerate(sources):
                        initLs = [bpy.data.objects[i.name]["Minimal Wavelength"] for i in sources] 
                        finalLs = [bpy.data.objects[i.name]["Maximal Wavelength"] for i in sources]
                    lengths = []
                    N = [bpy.data.objects[i.name]["Beams per ray"] for i in sources ]
           
                    for i in zip(N,zip(initLs,finalLs)):
                        
                        lengths.append(np.linspace(i[1][0],i[1][1],i[0]))
                    
                    for f in fsegs: 
                        for n,d in enumerate(lengths):
                            dirsMap = {}
                            for l in d:
                                dirsMap[l] = matrixData[1][n] 
                            dirs[str(n)+"_"+str(f)] = dirsMap
                    for n,i in enumerate(sources):
                        if bpy.data.scenes["Scene"].options.animating == True:
                            for frame in fsegs:
                                for src in lengths:
                                    for len in src:
                                        map = {}
                                        
                                        map[len] = 0
                                project[str(n)+"_"+str(frame)]=map          
                    originLoc = matrixData[0]
                    ###########################print(originLoc)
                    return dispersion(matrixData[1],targets,matrixData[0],S,aNs,nd,defs,sources,True,iors,frameselect,delta)
                    
            else:
                for n,i in enumerate(sources):
                    initLs = [bpy.data.objects[i.name]["Minimal Wavelength"] for i in sources] 
                    finalLs = [bpy.data.objects[i.name]["Maximal Wavelength"] for i in sources]
                lengths = []
                N = bpy.data.objects["Empty"]["Beams per ray"]
                for i in zip(N,zip(initLs,finalLs)):
                    lengths.append(np.linspace(i[1][0],i[1][1],i[0]))
                for n,i in enumerate(sources):
                    if bpy.data.scenes["Scene"].options.animating == True:
                        for wave in lengths:
                            map = {}
                            map[wave] = 0
                        project[str(n)]=map
                originLoc = matrixData[0]
                return dispersion(matrixData[1],targets,matrixData[0],S,aNs,nd,defs,sources,False,iors,bpy.context.scene.frame_current,delta)    
            
        
        
    else:
        
        originLoc = matrixData[0]
        ######################print(originLoc,matrixData[0])
        return dispersion(matrixData[1],targets,matrixData[0],S,aNs,nd,defs,sources,True,iors=iors,delta = delta)
        
       

        
##############MENU##############
class currIncCount(bpy.types.Operator):
    bl_idname = "currinccounter.change"
    bl_label = "alter the panel"
    
    def execute(self, context):
        try:
            bpy.context.scene["curveCount"]
        except:
            bpy.context.scene["curveCount"] = 0
            bpy.data.scenes['Scene']['curveDict'] = {}
            bpy.data.scenes['Scene']['current'] = {}

        for i in range(bpy.data.scenes["Scene"]["curveCount"]):
            setattr(bpy.types.Scene, f"chargecurve_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"current_{i}",bpy.props.FloatProperty(name = "Current"))
        bpy.data.scenes["Scene"]["curveCount"] += 1
        return {"FINISHED"}
    
class currDecCount(bpy.types.Operator):
    bl_idname = "currdeccounter.change"
    bl_label = "alter the panel"

    
    def execute(self, context):
        try:
            bpy.context.scene["curveCount"]
        except:
            bpy.context.scene["curveCount"] = 0
            bpy.data.scenes['Scene']['curveDict'] = {}
            bpy.data.scenes['Scene']['current'] = {}
        for i in range(bpy.data.scenes["Scene"]["curveCount"]):
            setattr(bpy.types.Scene, f"chargecurve_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"current_{i}",bpy.props.FloatProperty(name = "Current"))
        bpy.data.scenes["Scene"]["curveCount"] -= 1
        return {"FINISHED"}
class chargedIncCount(bpy.types.Operator):
    bl_idname = "chargeinccounter.change"
    bl_label = "alter the panel"
    
    def execute(self, context):
        try:
            bpy.context.scene["chargePrismCount"]
        except:
            bpy.context.scene["chargePrismCount"] = 0
            bpy.data.scenes['Scene']['chargeprismDict'] = {}
            bpy.data.scenes['Scene']['permittivity'] = {}

        for i in range(bpy.data.scenes["Scene"]["chargePrismCount"]):
            setattr(bpy.types.Scene, f"chargeprism_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"permit_{i}",bpy.props.FloatProperty(name = "Relative Electrical Permittivity", soft_min = 0.0001))
        bpy.data.scenes["Scene"]["ChargePrismCount"] += 1
        return {"FINISHED"}
    
class chargedDecCount(bpy.types.Operator):
    bl_idname = "chargedeccounter.change"
    bl_label = "alter the panel"

    
    def execute(self, context):
        try:
            bpy.context.scene["ChargePrismCount"]
        except:
            bpy.context.scene["ChargePrismCount"] = 0
            bpy.data.scenes['Scene']['ChargeprismDict'] = {}
            bpy.data.scenes['Scene']['permittivity'] = {}

        
        for i in range(bpy.data.scenes["Scene"]["MagPrismCount"]):
            setattr(bpy.types.Scene, f"chargeprism_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"permit_{i}",bpy.props.FloatProperty(name = "Relative Electrical Permittivity", soft_min = 0.0001))
        bpy.data.scenes["Scene"]["ChargePrismCount"] -= 1
        return {"FINISHED"}
    
class MagDecCount(bpy.types.Operator):
    bl_idname = "magdeccounter.change"
    bl_label = "alter the panel"
    
    def execute(self, context):
        try:
            bpy.context.scene["MagPrismCount"]
        except:
            bpy.context.scene["MagPrismCount"] = 0
            bpy.data.scenes['Scene']['MagprismDict'] = {}
            bpy.data.scenes['Scene']['permeability'] = {}

        for i in range(bpy.data.scenes["Scene"]["MagPrismCount"]):
            setattr(bpy.types.Scene, f"magprism_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"perm_{i}",bpy.props.FloatProperty(name = "Magnetic Permeability", soft_min = 0.0001))
        bpy.data.scenes["Scene"]["MagPrismCount"] -= 1
        return {"FINISHED"}
    
class MagIncCount(bpy.types.Operator):
    bl_idname = "maginccounter.change"
    bl_label = "alter the panel"

    
    def execute(self, context):
        try:
            bpy.context.scene["MagPrismCount"]
        except:
            bpy.context.scene["MagPrismCount"] = 0
            bpy.data.scenes['Scene']['MagprismDict'] = {}
            bpy.data.scenes['Scene']['permeability'] = {}

        
        for i in range(bpy.data.scenes["Scene"]["MagPrismCount"]):
            setattr(bpy.types.Scene, f"magprism_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"perm_{i}",bpy.props.FloatProperty(name = "Relative Magnetic Permeability", soft_min = 0.0001))
        bpy.data.scenes["Scene"]["MagPrismCount"] += 1
        return {"FINISHED"}
    
class DecCount(bpy.types.Operator):
    bl_idname = "deccounter.change"
    bl_label = "alter the panel"
    
    def execute(self, context):
        try:
            bpy.context.scene["PrismCount"]
        except:
            bpy.context.scene["PrismCount"] = 0
            bpy.data.scenes['Scene']['prismDict'] = {}
            bpy.data.scenes['Scene']['AbbesDict'] = {}

        for i in range(bpy.data.scenes["Scene"]["PrismCount"]):
            setattr(bpy.types.Scene, f"prism_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"aN_{i}",bpy.props.FloatProperty(name = "Abbe Number", soft_min = 0.0001))
            setattr(bpy.types.Scene, f"IOR_{i}",bpy.props.FloatProperty(name = "IOR", soft_min = 0))
        bpy.data.scenes["Scene"]["PrismCount"] -= 1

        return {"FINISHED"}
class IncCount(bpy.types.Operator):
    bl_idname = "inccounter.change"
    bl_label = "alter the panel"

    
    def execute(self, context):
        try:
            bpy.data.scenes['Scene']["PrismCount"]
        except:
            bpy.data.scenes['Scene']["PrismCount"] = 0
            bpy.data.scenes['Scene']['prismDict'] = {}
            bpy.data.scenes['Scene']['AbbesDict'] = {}

        
        for i in range(bpy.data.scenes["Scene"]["PrismCount"]):
            setattr(bpy.types.Scene, f"prism_{i}", bpy.props.PointerProperty(type=bpy.types.Object))
            setattr(bpy.types.Scene, f"aN_{i}",bpy.props.FloatProperty(name = "Abbe Number", soft_min = .1))
            setattr(bpy.types.Scene, f"IOR_{i}",bpy.props.FloatProperty(name = "IOR", soft_min = 0))
        bpy.data.scenes["Scene"]["PrismCount"] += 1
        return {"FINISHED"}

           
class ArrowEmpty(bpy.types.Operator):
    bl_idname = "add.empty_arrow"
    bl_label = "add empty arrow"

    def execute(self, context):  
        self.count = bpy.data.scenes["Scene"]["PrismCount"]
        bpy.ops.object.empty_add(type='SINGLE_ARROW', align='WORLD', location=(0, 0, 0), scale=(10, 10, 10))
        bpy.ops.transform.resize(value=(20, 20, 20), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', constraint_axis=(True, True, True), mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=False, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
        bpy.context.object["Type"] = "Ray Source"
        bpy.context.object["Beams per ray"] = 0
        
        ui = bpy.context.object.id_properties_ui("Beams per ray")
        ui.update(description = "Beams generated per ray")
        #ui.update(Type= "float")
        ui.update(default = 0)
        ui.update(min=0, soft_min=0)
        ui.update(max=200000000, soft_max=200000000)
        bpy.context.object["Minimal Wavelength"] = 0.4
        ui = bpy.context.object.id_properties_ui("Minimal Wavelength")
        ui.update(description = "Minimum wavelength of the ray bandwidth")
        ui.update(default = 0.4)
        ui.update(min=0.4, soft_min=0.4)
        ui.update(max=0.7, soft_max=0.7)
        bpy.context.object["Maximal Wavelength"] = 0.7
        ui = bpy.context.object.id_properties_ui("Maximal Wavelength")
        ui.update(description = "Maximum wavelength of the ray bandwidth")
        ui.update(default = 0.7)
        ui.update(min=0.4, soft_min=0.4)
        ui.update(max=0.7, soft_max=0.7)
        return {"FINISHED"}
    
class ElectricityMenu(bpy.types.Panel):
    bl_region_type = "UI"
    bl_idname="Electricity"
    bl_label="Electricitiy"
    bl_space_type="VIEW_3D"
    def draw(self,context):
        row = self.layout.row()
        row.operator("chargeinccounter.change",text="Add Charged Polygon")
        row.operator("chargedeccounter.change",text="Remove Charged Polyon")
        row.operator("currinccounter.change",text="Add Current")
        row.operator("currdeccounter.change",text="Remove Current")
        layout = self.layout
        bpy.data.scenes["Scene"]["ChargePrismCount"]
        for i in range(bpy.data.scenes["Scene"]["ChargePrismCount"]):
            try:
                self.layout.prop(context.scene, f'chargeprism_{i}')
                self.layout.prop(bpy.data.scenes["Scene"][f"chargeprism_{i}"],"rotation_euler")
                self.layout.prop(bpy.data.scenes["Scene"][f"chargeprism_{i}"],"location")
                self.layout.prop(context.scene,f'permit_{i}')

            except:
                pass    
        for i in range(bpy.data.scenes["Scene"]["curveCount"]):
            try:
                self.layout.prop(context.scene, f'currentCurve_{i}')
                self.layout.prop(context.scene,f'current_{i}')
            except:
                pass
class MagnetismMenu(bpy.types.Panel):
    bl_idname= "Magnetism"
    bl_label = "Magnetism"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"

    def draw(self,context):
        
        layout = self.layout
        self.layout.prop(context.scene.options,"Span")
        self.layout.prop(context.scene.options,"Resolution")
        row = self.layout.row()
        row.operator("maginccounter.change",text="Add Magnetic Prism")
        row.operator("magdeccounter.change",text="Remove Magnetic Prism")
        row.operator("timed.mag",text="Calculate Fields")
        
        bpy.data.scenes["Scene"]["MagPrismCount"]
        for i in range(bpy.data.scenes["Scene"]["MagPrismCount"]):
            try:
                self.layout.prop(context.scene, f'magprism_{i}')
                self.layout.prop(bpy.data.scenes["Scene"][f"magprism_{i}"],"rotation_euler")
                self.layout.prop(bpy.data.scenes["Scene"][f"magprism_{i}"],"location")
                self.layout.prop(context.scene,f'perm_{i}')
            except:
                pass
    
    
    
class RenderPanel(bpy.types.Panel):
    bl_idname = "RenderMenu"
    bl_label = "Rendering"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    
    def draw(self, context):
        layout = self.layout
        layout.prop(context.scene.options,"photonRadius")
        layout.prop(context.scene.options,"photonLength")
class AnimPanel(bpy.types.Panel):
    bl_idname = "AnimPanel"
    bl_label = "Animation"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    
    def draw(self, context):
        layout = self.layout
        self.layout.prop(context.scene.options,"animating")
        self.layout.prop(context.scene.options,"startTime")
        self.layout.prop(context.scene.options,"endTime")
        self.layout.prop(context.scene.options,"photonFreq")
        self.layout.prop(context.scene.options,"speed")
class VIEW3D_PT_my_custom_panel(bpy.types.Panel):
    pass
    bl_space_type="VIEW_3D"
    bl_region_type="UI"
    bl_label="Light"
    def draw(self,context):
        row = self.layout.row()
        row.operator("timed.ray",text="Calculate Rays")
        row.operator("add.empty_arrow",text="Create Source")
        row.operator("inccounter.change",text="Add prism")
        row.operator("deccounter.change",text="Remove prism")
        self.layout.prop(context.scene.options,"RecursiveDepth")
        self.layout.menu("AnimPanel")
        self.layout.menu("RenderMenu")
        scene = context.scene
        
        for i in bpy.data.objects:
            
            try:
                if i["Type"] == "Ray Source":
                    self.layout.prop(i,"name")
                    self.layout.prop(i,"location")
                    self.layout.prop(i,"rotation_euler")
                    self.layout.prop(i,'["Maximal Wavelength"]')
                    self.layout.prop(i,'["Minimal Wavelength"]')
                    self.layout.prop(i,'["Beams per ray"]')
            except:
                pass
     
        try:
            bpy.data.scenes["Scene"]["PrismCount"]
            for i in range(bpy.data.scenes["Scene"]["PrismCount"]):
                self.layout.prop(context.scene, f'prism_{i}')
                self.layout.prop(context.scene,f'aN_{i}')
                self.layout.prop(context.scene,f'IOR_{i}')
        except:
            pass
class options(bpy.types.PropertyGroup):
    Resolution: bpy.props.IntProperty(name = "Domain Resolution",min = 0)
    Span: bpy.props.IntProperty(name ="Field Span", soft_min = 0)
    RecursiveDepth: bpy.props.IntProperty(name = "Depth", soft_min=0)
    startTime: bpy.props.IntProperty(name="Initial Frame")
    endTime: bpy.props.IntProperty(name = "Final Frame")
    animating: bpy.props.BoolProperty(name = "animating")
    photonFreq: bpy.props.IntProperty(name = "Photon Frequency", min=0)
    speed: bpy.props.FloatProperty(name = "Photon Speed", default=299792458, min = 0)
    photonLength: bpy.props.FloatProperty(name = "Photon Length", default=1, min = 0)
    photonRadius: bpy.props.FloatProperty(name = "Photon Cross Section Radius", default=1, min = 0)
classes = [options,RenderPanel,AnimPanel,ArrowEmpty,IncCount,currDecCount,currIncCount,DecCount,MagDecCount,MagIncCount,MagnetismMenu,chargedDecCount,chargedIncCount,ElectricityMenu,MagModal,VIEW3D_PT_my_custom_panel,Modal]

def register():
    for cls in classes:
        bpy.utils.register_class(cls)
        bpy.utils.register_class(cls)
        bpy.types.Scene.options = bpy.props.PointerProperty(type = options)
def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)
        del bpy.types.Scene.prop
        del bpy.types.Object.Arrow
        del bpy.types.Scene.options
if __name__ =="__main__":
    register()


