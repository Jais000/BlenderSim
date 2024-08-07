import bpy



class Currents:
    def __init__(self):
        self.currents = [bpy.data.scenes["Scene"][f"currentCurve_{i}"] for i in range(bpy.data.scenes["Scene"]["curveCount"]) if bpy.data.scenes["Scene"][f"currentCurve_{i}"] and bpy.data.scenes["Scene"][f"currentCurve_{i}"] in [l for l in bpy.data.scenes["Scene"].objects]] 
        self.dirDict = {}
    def current(self):   
        try:
            bpy.context.view_layer.objects.active.select_set(False)
        except:
            pass
        for i in self.currents:
            self.dirDict[i.name] = {}
            mesh = bpy.data.meshes.new_from_object(i)
            new_obj = bpy.data.objects.new(f"{i.name}_copy", mesh)
            new_obj.matrix_world = i.matrix_world
            bpy.context.collection.objects.link(new_obj)
            verts = [(new_obj.matrix_world @ v.co) for v in new_obj.data.vertices]
            for j in range(len(verts)):
                if j < len(verts)-1:
                    self.dirDict[i.name][f"{verts[j].x},{verts[j].y},{verts[j].z}"] = verts[j+1]- verts[j]
            new_obj.select_set(False)
        #print(self.dirDict)
    def clean(self):
        data = []
        for i in bpy.context.scene.objects:
            objs = bpy.data.objects
            if "_copy" in i.name:
                data.append(i.data)
                objs.remove(i, do_unlink=True)
        for i in data:
            bpy.data.meshes.remove(i) 


def chargedPolies():
    polies = [bpy.data.scenes["Scene"][f"chargeprism_{i}"] for i in range(bpy.data.scenes["Scene"]["chargePrismCount"]) if bpy.data.scenes["Scene"][f"chargeprism_{i}"] and bpy.data.scenes["Scene"][f"chargeprism_{i}"] in [l for l in bpy.data.scenes["Scene"].objects]] 

