import bpy, bmesh
import numpy as np
from mathutils import Vector 
import pyopenvdb as vdb
from numpy import sqrt
from mathutils import Euler
import os
import sys
dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir)
from Electricity import Currents
import math 
val = {}
vel = {}
force = {}
angularVel = {}
#bpy.ops.mesh.primitive_cylinder_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
obj = bpy.context.object

def magForce(r,m1,m2):
    mm1 = Vector((0,0,1))
    mm1.rotate(m1.rotation_euler)
    mm2 = Vector((0,0,1))
    mm2.rotate(m2.rotation_euler)
    mu = 1.25663706* 10**(-6)*10**4
    absr = sqrt(r[0]**2+r[1]**2+r[2]**2)
    a=((3*mu)/(4*np.pi*(absr**5)))
    b1 = (mm1.dot(r))*mm2
    b2 = (mm2.dot(r))*mm1
    b3 = (mm1.dot(mm2))*r
    b4 = Vector(((  (5*(mm1.dot(r))*(mm2.dot(r)) )/( absr**2 ))*r))
    force = Vector((a*(b1+b2+b3-b4)))
    print("ACC",force)
    return force
def magTorque(mass,mu,m2,m1,r):
    mm1 = Vector((0,0,1))
    mm1.rotate(m1.rotation_euler)
    mm2 = Vector((0,0,1))
    mm2.rotate(m2.rotation_euler)    
    absr = sqrt(r[0]**2+r[1]**2+r[2]**2)
    l = 10
    mp = 1.25663706* 10**(-6)*10**6
    a = mp/(4*np.pi)
    b = (3*r.normalized()*(r.normalized().dot(mm1)))-mm1
    c = absr**3
    B = a*(b/c)
    t = mm2.cross(B)
    acc = t/(mass*l**2)
    return acc


def clear():
    for i in bpy.context.scene.objects:
            objs = bpy.data.objects
            if "crossSection" in i.name:  
                objs.remove(i, do_unlink=True)
                
    for i in bpy.data.meshes:
            objs = bpy.data.objects
            if "_crossSection_Mesh" in i.name:
                bpy.data.meshes.remove(i) 
def getMagneticArea(obj):
    new_obj = obj.copy()
    new_obj.data = obj.data.copy()
    new_obj.data.name =  f"{obj.name}_crossSection_Mesh"
    new_obj.name = f"{obj.name}_crossSection"
    new_obj.matrix_world = obj.matrix_world
    loc = new_obj.location
    rot = new_obj.matrix_world.to_euler()
    bpy.context.scene.collection.objects.link(new_obj)
    bpy.ops.mesh.primitive_plane_add(location=new_obj.location, scale=(100, 100, 100))
    bpy.context.object.name="slicer" 
    objVerts = [(obj.matrix_world @ v.co) for v in obj.data.vertices]
    start = 0
    num = 0 
    for i in objVerts:
        dis = sqrt((i[0]-loc[0])**2+(i[1]-loc[1])**2 + (i[2]-loc[2])**2)
        start = max(dis,start)
        if start == dis:
            ######print("TRYE",dis)
            past = max(max(abs(i[0]),abs(i[1])),abs(i[2]))
            num = max(num,abs(past))
    bpy.ops.transform.resize(value=(num, num, num), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    bpy.context.scene.objects['slicer'].rotation_euler = new_obj.matrix_world.to_euler() 
    mod = new_obj.modifiers.new(type="BOOLEAN", name="cross")
    mod.object = bpy.context.scene.objects['slicer']
    mod.operation = "INTERSECT"
    new_obj.select_set(True)
    bpy.context.view_layer.objects.active = new_obj
    bpy.ops.object.modifier_apply(modifier="cross")   
    for i in bpy.context.scene.objects:
        objs = bpy.data.objects
        if "slicer" in i.name:  
            data = i.data
            objs.remove(i, do_unlink=True)
            bpy.data.meshes.remove(data)
    return new_obj 
    
def createFieldPoint(rfd,dir,r,prism):
        
        try:
            n = prism.split("_")[1]
            mu = 1.25663706* 10**(-6)*bpy.data.scenes["Scene"][f"perm_{n}"]
        except:
            mu = 1.25663706* 10**(-6)
       
        volume = targets[prism]['volume']
        area = targets[prism]['mArea']  
        m = (rfd*volume)/(1.25663706* 10**(-6))
        i = m/area
        print("TST",volume,area,m,i)
        a= (i*mu)/(4*np.pi)
        b= (dir.cross(r.normalized()))/((sqrt(r[0]**2+r[1]**2+r[2]**2))**3)
        print("CC",a*b)
        return a*b
    
def deleteField():
    for i in bpy.context.scene.objects:
            objs = bpy.data.objects
            if "MagFieldVec" in i.name:  
                objs.remove(i, do_unlink=True)
                continue
            if "domain" in i.name:  
                objs.remove(i, do_unlink=True)
                continue
            if "Force" in i.name:
                objs.remove(i,do_unlink=True)
                
    for i in bpy.data.meshes:
            objs = bpy.data.objects
            if "Grid" in i.name:
                bpy.data.meshes.remove(i) 
                continue
            if "Circle" in i.name:
                bpy.data.meshes.remove(i) 
              

def createGrid(obj,resolution):
    start = 0
    num = 0 
    past = 0
    loc = obj.location
    objVerts = [(obj.matrix_world @ v.co) for v in obj.data.vertices]
    for i in objVerts:
        dis = sqrt((i[0]-loc[0])**2+(i[1]-loc[1])**2 + (i[2]-loc[2])**2)
        start = max(dis,start)
        if start == dis:
            ######print("TRYE",dis)
            past = max(max(abs(i[0]),abs(i[1])),abs(i[2]))
            num = max(num,abs(past))
    ###print("NUM",num)
    domain = np.arange(-num,num,resolution)
    ######print("DOMAIN",objVerts)
    for i in range(len(domain)):
        bpy.ops.mesh.primitive_grid_add(enter_editmode=True, align='WORLD', location=(loc[0],loc[1],loc[2]-domain[i]), scale=(4, 4, 4))
        bpy.ops.transform.resize(value=(num, num, num), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
    #bpy.context.view_layer.objects.active = bpy.context.object
    bpy.ops.object.mode_set(mode="OBJECT")
    mod = bpy.context.object.modifiers.new(type="BOOLEAN", name="confine")
    mod.object = obj
    mod.operation = "INTERSECT"
    bpy.ops.object.modifier_apply(modifier="confine")
    print("CREATEGRID",obj)
    return [bpy.context.object,obj]
    
def createMagnet(args,attract):
    
    gridVerts = [(args[0].matrix_world @ v.co) for v in args[0].data.vertices]
    rads = Vector((0,0,0))
    if attract == True:
        rotation = args[1].rotation_euler
        bpy.ops.object.empty_add(type='SINGLE_ARROW', align='WORLD', location=args[1].location, scale=(1, 1, 1))
        bpy.context.object.name = f"domain_{args[1].name}"
        bpy.context.object.rotation_euler = rotation
        args[0].hide_set(True)
        args[0].select_set(False)
        args[1].select_set(True)
        bpy.context.view_layer.objects.active = args[1]

        #bpy.context.scene.objects.active.rigid_body.type=True
        #bpy.data.objects["Cube.001"].field.type
    else:
        rotation = args[1].rotation_euler  
        s = bpy.data.scenes["Scene"]["options"]["Span"]
        createPoints(rotation,[1,1,1],(s,s,s),args[1])
        args[0].hide_set(True)

def createPoints(dir,m,span,object):
    for i in [k for k in bpy.data.scenes["Scene"].keys()]:
        if "magprism" in i:
            if bpy.data.scenes["Scene"][f"{i}"] == object:
                prism = i
    global val
    val = {} 
    vectorSource = {}
    #dir = Vector((1,1,1))
    #diff = dir.rotation_difference((0,0,1)).to_euler("XYZ")
    bpy.ops.mesh.primitive_circle_add(enter_editmode=False, align='WORLD', location=object.location, scale=(1, 1, 1))
    obj = bpy.context.object
    obj.name = f"MagCircle_{object.name}"
    obj.rotation_euler = dir
    bpy.context.scene.objects[f"MagCircle_{object.name}"].select_set(False)
    bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
    objVerts = [(obj.matrix_world @ v.co) for v in obj.data.vertices]

    for i in range(0,len(objVerts)):
        if not i >= len(objVerts)-1:
            vectorSource[i] = {'loc':objVerts[i],'dir':objVerts[i]-objVerts[i+1]}

    for i in range(-span[0],span[0]):
       for j in range(-span[1],span[1]):
           for k in range(-span[2],span[2]):
                print(i,j,k,object)
                for data in vectorSource.values():
                    delr =Vector((i,j,k))-data['loc']
                    try:
                        try: 
                            val[f'{i},{j},{k}'] += createFieldPoint(1,data['dir'],delr,object) 
                        except:
                            val[f'{i},{j},{k}'] = createFieldPoint(1,data['dir'],delr,object) 
                    except:
                        try:
                            val[f'{i},{j},{k}'] += Vector((0,0,0))
                        except:
                            val[f'{i},{j},{k}'] = Vector((0,0,0))
                
def createCurrentPoints():
    span = bpy.context.scene["options"]["Span"]
    electricity = Currents()
    electricity.clean()
    electricity.current()    

    for i in range(-span,span):
       for j in range(-span,span):
           for k in range(-span,span):
                for n,l in enumerate(electricity.dirDict):
                    for data in electricity.dirDict[l]:
                        
                        vec = Vector((float(data.split(',')[0]),float(data.split(',')[1]),float(data.split(',')[2])))
                        try:
                            print(bpy.data.scenes["Scene"][f"current_{n}"])
                            print("SUCCESS") 
                            try: 
                                val[f'{i},{j},{k}'] += createFieldPoint(bpy.data.scenes["Scene"][f"current_{n}"],electricity.dirDict[l][data],Vector((i,j,k))-vec,False) 

                            except:
                                val[f'{i},{j},{k}'] = createFieldPoint(bpy.data.scenes["Scene"][f"current_{n}"],electricity.dirDict[l][data],Vector((i,j,k))-vec,False) 
                        except:
                            try:
                                val[f'{i},{j},{k}'] += Vector((0,0,0))
                                print("TEST2")
                            except:
                                val[f'{i},{j},{k}'] = Vector((0,0,0))
                                print("TEST4")                    
def createForce():
    magnets = []
    for i in bpy.context.scene.objects:
        if "domain" in i.name:
            name= f"{i.name.split('_')[1].split('.')[0]}"
            if not name in [i.name for i in magnets]:
                magnets.append(bpy.context.scene.objects[name])
    for n,k in enumerate(magnets):
        mass = 1
        force[k] = {}
        for i in [j for j in bpy.context.scene.objects if f"domain_{k.name}" in j.name]:
            force[k]['linForce'] = Vector((0,0,0))
            force[k]['rotAcc'] = Vector((0,0,0))
            acc = Vector((0,0,0))
            for p in [o for o in bpy.context.scene.objects if f"domain" in o.name and not k.name in o.name]:

                force[k]['linForce'] += magForce(i.location-p.location,i,p)
                force[k]['rotAcc'] += magTorque(mass, bpy.data.scenes["Scene"][f"perm_{n}"],p,i,i.location-p.location)
            #s = sqrt(f[0]**2+f[1]**2+f[2]**2)

            #bpy.ops.transform.resize(value=(s,s,s), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
def animate(frame):
    global angularVel
    global vel
    magnets = []
    for i in bpy.context.scene.objects:
        if "domain" in i.name:
            name= f"{i.name.split('_')[1].split('.')[0]}"
            if not name in [i.name for i in magnets]:
                magnets.append(bpy.context.scene.objects[name])
    mass = 1
    for k in magnets:
        
        magnet = k
        bpy.context.view_layer.objects.active = magnet
        if frame == bpy.context.scene.frame_start:
            bpy.ops.rigidbody.object_add()
            magnet.rigid_body.collision_shape = 'SPHERE'
            magnet.rigid_body.friction = 0.6
            magnet.rigid_body.restitution = 0
            magnet.rigid_body.use_margin = True
            magnet.rigid_body.kinematic=True
        fps = bpy.context.scene.render.fps
        #print(fps,frame)
        time = frame/fps
        
        
        
        #pref_edit = bpy.context.user_preferences.edit
        #keyInterp = pref_edit.keyframe_new_interpolation_type
        #pref_edit.keyframe_new_interpolation_type ='LINEAR'
        start_position = magnet.location
        ##print("START",magnet.location)
        
        acc = (force[magnet]['linForce']/mass)
        print(acc,magnet.name)
        rotacc = force[magnet]['rotAcc']
        try:
            angularVel[magnet] += rotacc
        except:
            angularVel[magnet] = rotacc
        try:
            vel[magnet] += acc
        except:
            vel[magnet] = acc
        omegaVel = Euler(angularVel[magnet], "XYZ")
        o = magnet.rotation_euler 
        bpy.context.scene.objects[f"domain_{magnet.name}"].rotation_euler = Euler((o.x + omegaVel.x,o.y + omegaVel.y,o.z + omegaVel.z),"XYZ")
        bpy.context.scene.objects[f"domain_{magnet.name}"].location = start_position + vel[magnet]
        magnet.rotation_euler = Euler((o.x + omegaVel.x,o.y + omegaVel.y,o.z + omegaVel.z),"XYZ")
        magnet.location = start_position + vel[magnet]
        magnet.keyframe_insert(data_path='location', frame=frame)
        magnet.keyframe_insert(data_path='rotation_euler',frame=frame)
        bpy.context.scene.objects[f"domain_{magnet.name}"].keyframe_insert(data_path='location', frame=frame)
        bpy.context.scene.objects[f"domain_{magnet.name}"].keyframe_insert(data_path='rotation_euler',frame=frame)
        ##print("END",magnet.location)
        
        #pref_edit.keyframe_new_interpolation_type = keyInterp
angV = Euler((0,0,0), "XYZ")
def createField():
    s = bpy.data.scenes["Scene"]["options"]["Span"]
    span = [s,s,s]
    global val

    for i in range(-span[0],span[0]):
       for j in range(-span[1],span[1]):
           for k in range(-span[2],span[2]):
                
                dist = val[f'{i},{j},{k}'].magnitude
                print("TEST",dist)
                try:
                    dist = math.log(dist*1000,20)
                except:
                    dist = 0
                val[f'{i},{j},{k}'] = val[f'{i},{j},{k}'].normalized()
                bpy.ops.object.empty_add(type='SINGLE_ARROW', align='WORLD', location=(i,j,k), scale=(1,1,1))
                bpy.ops.transform.resize(value=(dist, dist, dist), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=False, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, snap=False, snap_elements={'INCREMENT'}, use_snap_project=False, snap_target='CLOSEST', use_snap_self=True, use_snap_edit=True, use_snap_nonedit=True, use_snap_selectable=False)
                arr = bpy.context.object
                arr.name = "MagFieldVec"
                eul = val[f'{i},{j},{k}'].to_track_quat("-Z","Z").to_euler("XYZ")
                arr.rotation_euler = eul
                arr.select_set(False)
deleteField()
clear()
print("VAL",val)       
targets ={}
for i in range(bpy.data.scenes["Scene"]["MagPrismCount"]) :
    if bpy.data.scenes["Scene"][f"magprism_{i}"] and bpy.data.scenes["Scene"][f"magprism_{i}"] in [l for l in bpy.data.scenes["Scene"].objects]:
        targets[bpy.data.scenes["Scene"][f"magprism_{i}"]] = {}
        bm = bmesh.new()
        bm.from_mesh( bpy.data.scenes["Scene"][f"magprism_{i}"].data )
        volume = float(bm.calc_volume())
        targets[bpy.data.scenes["Scene"][f"magprism_{i}"]]["volume"]=volume
        crossSection = getMagneticArea(bpy.data.scenes["Scene"][f"magprism_{i}"])
        bm = bmesh.new()
        bm.from_mesh( crossSection.data )
        area = sum(f.calc_area() for f in bm.faces)
        targets[bpy.data.scenes["Scene"][f"magprism_{i}"]]["mArea"]=area
print(targets)

for i in targets:
    createMagnet(createGrid(i,3),True)
#createCurrentPoints()
#createField()
    #for k in np.arange(0,300,1):
        #createForce()
        #animate(k)
#

class MagModal(bpy.types.Operator):
    bl_idname = "timed.mag"
    bl_label = "mag"
    
    def modal(self, context, event):
        global val    
        if event.type in {'RIGHTMOUSE', 'ESC'}:
            self.cancel(context)
            return {'CANCELLED'}
        
        cont = bpy.context.object
        if event.type == 'TIMER':
            res = bpy.data.scenes["Scene"].options.Resolution
            val = {}
            deleteField()
            targets =[bpy.data.scenes["Scene"][f"magprism_{i}"] for i in range(bpy.data.scenes["Scene"]["MagPrismCount"]) if bpy.data.scenes["Scene"][f"magprism_{i}"] and bpy.data.scenes["Scene"][f"magprism_{i}"] in [l for l in bpy.data.scenes["Scene"].objects]]
            for i in targets:
                createMagnet(createGrid(i,res),True)
            createCurrentPoints()
            createField()


        return {'PASS_THROUGH'}
    
    def execute(self, context):
        global val
        val = {}
        T = .1 
        wm = context.window_manager
        self._timer = wm.event_timer_add(T, window=context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    def cancel(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer) 
        

