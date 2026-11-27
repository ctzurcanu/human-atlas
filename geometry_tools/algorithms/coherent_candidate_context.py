"""Use one complete candidate snapshot for both sides of a contact check."""
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha


def context_snapshots(current,part_id,rows):
    row=rows.get(part_id)
    if row is None:return current,current
    if sha(row['path'])!=row['sha256']or sha(row['sourcePath'])!=row['sourceSha256']:
        raise ValueError('Changed context must bind its completed source and candidate files')
    source,new=Mesh.load(row['sourcePath']),Mesh.load(row['path'])
    if source.geometry_sha!=current.geometry_sha or source.attributes_sha!=current.attributes_sha:
        raise ValueError('Changed context source is stale relative to the current catalogue')
    if new.frame!=current.frame or not np.array_equal(new.faces,current.faces):
        raise ValueError('This source-preserving contact stage requires unchanged peer face ownership')
    return current,new
