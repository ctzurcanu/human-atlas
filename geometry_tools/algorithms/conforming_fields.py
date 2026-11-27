"""Replay source-corner fields on an explicitly owned Boolean arrangement."""
import numpy as np


def replay_corner_fields(positions, surfaces, normals, uvs, result):
    q = result['positions']
    output_normals, output_uvs = {}, {}
    for sid, fs in result['surfaces'].items():
        chart = result['charts'][sid]
        refs = np.asarray(chart['source_surfaces'])
        parents, signs, bary = chart['source_faces'], chart['normal_signs'], chart['barycentrics']
        words = np.empty((len(fs), 3, 3), '<i2')
        tex = np.empty((len(fs), 3, 2), '<f4')
        for family in np.unique(refs):
            ids = np.flatnonzero(refs == family)
            par, w = parents[ids], bary[ids]
            original = normals[family][par]
            values = np.einsum('tik,tkj->tij', w, original.astype(float))
            length = np.linalg.norm(values, axis=2)
            if np.any(length == 0):
                raise ValueError('A source normal chart interpolates to zero')
            ww = np.rint(values / length[:, :, None] * 32767).astype('<i2')
            tt = np.einsum('tik,tkj->tij', w, uvs[family][par]).astype('<f4')
            exact = np.all(q[fs[ids]][:, :, None] == positions[surfaces[family][par]][:, None], axis=3)
            for corner in range(3):
                rows = np.flatnonzero(exact[:, corner].any(1))
                match = exact[:, corner].argmax(1)
                ww[rows, corner] = original[rows, match[rows]]
                tt[rows, corner] = uvs[family][par[rows], match[rows]]
            words[ids] = ww * signs[ids, None, None]
            tex[ids] = tt
        output_normals[sid], output_uvs[sid] = words, tex
    return output_normals, output_uvs
