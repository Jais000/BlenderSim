import bpy
import sys
import os

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir)
import refraction as light

class VIEW3D_PT_my_custom_panel(bpy.types.Panel):
    pass
    bl_space_type="VIEW_3D"
    bl_region_type="UI"
    bl_label="Light"
    def draw(self,context):
        row = self.layout.row()
        row.operator("timed.ray",text="Calculate Rays")
        scene = context.scene
        mytool = scene.my_tool
        self.layout.prop(mytool,"BandDefinition")
        self.layout.prop(mytool,"a_n")
        self.layout.prop(mytool,"Length_i")
        self.layout.prop(mytool,"Length_f")

        self.layout.prop(context.scene, "Prism_Collection")    
class options(bpy.types.PropertyGroup):
    BandDefinition: bpy.props.IntProperty(name="Rays Per Band",soft_min=0, soft_max = 1000)
    a_n: bpy.props.FloatProperty(name="Abbes Number",soft_min=0, soft_max = 1000)
    Length_i: bpy.props.FloatProperty(name="Initial Wavelength",soft_min=0, soft_max = 1000) 
    Length_f: bpy.props.FloatProperty(name="Final Wavelength",soft_min=0, soft_max = 1000) 
    

classes = [VIEW3D_PT_my_custom_panel,options]
def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.my_tool=bpy.props.PointerProperty(type=options)
    bpy.types.Scene.Prism_Collection = bpy.props.PointerProperty(type=bpy.types.Collection)


def unregister():
    light.unregister()
    for cls in classes:
        bpy.utils.unregister_class(cls)
    del bpy.types.Scene.my_tool
    del bpy.types.Scene.prop
if __name__ =="__main__":
    #light.register()
    register()

