import json
import numpy as np
import bpy

from mathutils import *

def vec_from_array(arr):
    return Vector((
        float(arr[0]),
        float(arr[1]),
        float(arr[2])
    ))

def arr_to_buffer(buffer, arr, i):
    buffer[i*3]   = float(arr[i][0])
    buffer[i*3+1] = float(arr[i][1])
    buffer[i*3+2] = float(arr[i][2])


def main():
    #print(os.getcwd())
    with open("/home/nick/Documents/School/Optimization/outputs/outp.json") as f:
        struct = json.load(f)
        
    N = len(struct['normals'])
    verts = []
    positions = np.zeros(3*N)
    directions = np.zeros(3*N)
    F = np.zeros(3*N)
    for i in range(N):
        verts.append(Vector((0, 0, 0)))
        arr_to_buffer(positions, struct['positions'], i)
        arr_to_buffer(directions, struct['normals'], i)
        arr_to_buffer(F, struct['F'], i)
        
    obj_name = "ACS"
    if obj_name in bpy.data.objects and obj_name+"_data" in bpy.data.meshes:
        # Object already exists
        obj = bpy.data.objects[obj_name]
        mesh_data = bpy.data.meshes[obj_name + "_data"]
        
    else:
        # Create new object
        mesh_data = bpy.data.meshes.new(obj_name + "_data")
        obj = bpy.data.objects.new(obj_name, mesh_data)
        bpy.context.scene.collection.objects.link(obj)
        mesh_data.from_pydata(verts, [], [])
        #obj.data.attributes.new(name='positions', type='FLOAT_VECTOR', domain='POINT')
        obj.data.attributes.new(name='directions', type='FLOAT_VECTOR', domain='POINT')
        obj.data.attributes.new(name='forces', type='FLOAT_VECTOR', domain='POINT')
        
    obj.data.attributes['position'].data.foreach_set('vector', positions)
    obj.data.attributes['directions'].data.foreach_set('vector', directions)
    obj.data.attributes['forces'].data.foreach_set('vector', F)
    

if __name__ == "__main__":
    main()