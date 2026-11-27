"""Transport source normal fields independently on explicitly owned sheets.

Coincident coordinates alone do not authorize pooling deformation stars across
distinct source sheets. Every face and used raw normal control needs exactly
one declared owner; shared raw controls require a different, explicit assembly.
"""
import numpy as np
from .surface_normals import transport_vertex_words


def transport_partitioned_vertex_words(before,after,faces,words,partitions,*,maximum_facet_rotation_degrees):
    before,after=np.asarray(before),np.asarray(after);faces=np.asarray(faces);words=np.asarray(words)
    assert before.shape==after.shape and words.shape==before.shape
    coverage=np.concatenate([np.asarray(p,int) for p in partitions])
    np.testing.assert_array_equal(np.sort(coverage),np.arange(len(faces)))
    owners=np.full(len(before),-1,int);result=words.copy();proofs=[]
    for number,ids in enumerate(partitions):
        ids=np.asarray(ids,int);used,local=np.unique(faces[ids],return_inverse=True)
        assert np.all(owners[used]<0),'Raw normal control belongs to multiple declared sheets'
        owners[used]=number
        transported,proof=transport_vertex_words(before[used],after[used],local.reshape(-1,3),words[used],
            maximum_facet_rotation_degrees=maximum_facet_rotation_degrees)
        result[used]=transported
        stationary=np.array_equal(before[used],after[used])
        if stationary:np.testing.assert_array_equal(result[used],words[used])
        proofs.append(dict(partition=number,sourceFaceIds=ids.tolist(),sourceRawControls=used.tolist(),
            allPartitionCoordinatesFixed=stationary,transport=proof))
    return result,dict(complete=True,explicitFaceOwnership=True,noCrossSheetNormalStarPooling=True,
        partitions=proofs,unusedNormalWordsLiteral=True,maximumFacetRotation_degrees=maximum_facet_rotation_degrees)
